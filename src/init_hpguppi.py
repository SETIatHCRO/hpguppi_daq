#!/usr/bin/env python3
import argparse
import os, sys
import yaml
import subprocess
import re
import socket
import datetime

default_prefix_exec =	"/opt/mnt/bin/"
default_prefix_lib 	=	"/opt/mnt/lib/"
default_logdir			= None # Switches off logging
############ Hands off from here on
def index_yaml_value(yaml_value, instance, cpu_count, subsystem):
	if isinstance(yaml_value, list):
		return yaml_value[instance]
	if isinstance(yaml_value, dict):
		if cpu_count in yaml_value:
			return yaml_value[cpu_count]
		return yaml_value[subsystem]
	return yaml_value


def collect_yaml_value(yaml_value, instance, cpu_count, subsystem):
	collection = []
	for value in yaml_value:
		contextual_value = index_yaml_value(value, instance, cpu_count, subsystem)
		if isinstance(contextual_value, list):
			# handle nested containers
			collection += list(
				map(
					lambda v: index_yaml_value(v, instance, cpu_count, subsystem) if isinstance(v, list) or isinstance(v, dict) else v,
					contextual_value
				)
			)
		else:
			collection.append(contextual_value)
	return collection


def get_cpu_core_count():
	# Gather cpu core count
	cores_per_cpu = subprocess.run('grep cpu.cores /proc/cpuinfo'.split(' '), capture_output=True).stdout

	cores_per_cpu = cores_per_cpu.decode().strip()
	m = re.match(r'cpu cores\s*:\s*(\d+).*', cores_per_cpu)
	return int(m.group(1))


def get_extended_environment(environment_keys, keyword_variable_dict):
	env = os.environ
	for env_kv in environment_keys:
		key, val = env_kv.split('=')

		if '$'+key in val:
			replacement = env[key] if key in env else ''
			val = val.replace('$'+key, replacement)
		for var,keyword_val in keyword_variable_dict.items():
				val = val.replace('${}'.format(var), str(keyword_val))

		env[key] = val
	return env


def delete(
	system_name,
	instance,
	yaml_config,
	environment_keys = [],
	dry_run = False,
	config_filename = '???'
):

	subsystem_split_index = system_name.find('_')

	system_full_name = system_name
	# Select configuration for the system
	if subsystem_split_index > -1 and system_name not in yaml_config:
		subsystem = system_name[subsystem_split_index+1:]
		system_name = system_name[0:subsystem_split_index]
		assert system_name in yaml_config , 'Root-system {} not defined in {}'.format(system_name, config_filename)
		print('Accessing root-system \'{}\''.format(system_name))
	else:
		subsystem = ''
		assert system_name in yaml_config , 'System {} not defined in {}'.format(system_name, config_filename)

	system = yaml_config[system_name]
	
	prefix_exec = system.get('prefix_exec', default_prefix_exec)

	# Gather cpu core count
	try:
		cores_per_cpu = get_cpu_core_count()
	except:
		print('Error trying to get the cpu core count.')
		exit(0)
	
	system_environment_keys_start_index = len(environment_keys)
	if 'hashpipe_keyfile' in system:
		environment_keys.append('HASHPIPE_KEYFILE={}'.format(system['hashpipe_keyfile']))
	if 'environment' in system:
		environment_keys += collect_yaml_value(system['environment'], instance, cores_per_cpu, subsystem)
	
	cmd = [
		f"{prefix_exec}hashpipe_clean_shmem",
		"-I", str(instance),
		"-d"
	]

	print(' '.join(cmd))
	hashpipe_env = get_extended_environment(environment_keys, {})
	print(subprocess.run(cmd, env=hashpipe_env, check=True, stdout=subprocess.PIPE).stdout.decode())

def run(
	system_name,
	instance,
	yaml_config,
	additional_arguments = [],
	options = [],
	environment_keys = [],
	dry_run = False,
	config_filename = '???'
):
	
	subsystem_split_index = system_name.find('_')

	system_full_name = system_name
	# Select configuration for the system
	if subsystem_split_index > -1 and system_name not in yaml_config:
		subsystem = system_name[subsystem_split_index+1:]
		system_name = system_name[0:subsystem_split_index]
		assert system_name in yaml_config , 'Root-system {} not defined in {}'.format(system_name, config_filename)
		print('Accessing root-system \'{}\''.format(system_name))
	else:
		subsystem = ''
		assert system_name in yaml_config , 'System {} not defined in {}'.format(system_name, config_filename)

	system = yaml_config[system_name]

	# Gather cpu core count
	try:
		cores_per_cpu = get_cpu_core_count()
	except:
		print('Error trying to get the cpu core count.')
		exit(0)

	# Gather system configuration for the cores_per_cpu
	assert 'cpu_core_count_config' in system, '{} not defined for system {} in {}'.format('cpu_core_count_config', system_name, config_filename)
	if cores_per_cpu not in system['cpu_core_count_config']:
		print('{}[{}] not defined for system {} in {}'.format('cpu_core_count_config', cores_per_cpu, system_name, config_filename))
		exit(1)
	system_ccc_config = system['cpu_core_count_config'][cores_per_cpu] if isinstance(system['cpu_core_count_config'], dict) else True

	# Gather numanode_bind specifications
	instance_numanode_bind = instance
	if 'instance_numanode_bind' in system:
		if system['instance_numanode_bind'] is False:
			instance_numanode_bind = False
		else:
			instance_numanode_bind = system['instance_numanode_bind']
			if isinstance(instance_numanode_bind, dict):
				instance_numanode_bind = instance_numanode_bind[cores_per_cpu]
			if isinstance(instance_numanode_bind, list):
				instance_numanode_bind = instance_numanode_bind[instance]
	else:
		print('{} not found for system {} in {}, numactl binding matches instance enumeration'.format('instance_numanode_bind', system_name, config_filename))

	instance_numanode_cpubind = instance_numanode_bind
	if 'instance_numanode_cpubind' in system:
		instance_numanode_cpubind = system['instance_numanode_cpubind']
		if isinstance(instance_numanode_cpubind, dict):
			instance_numanode_cpubind = instance_numanode_cpubind[cores_per_cpu]
		if isinstance(instance_numanode_cpubind, list):
			instance_numanode_cpubind = instance_numanode_cpubind[instance]

	instance_numanode_membind = instance_numanode_bind
	if 'instance_numanode_membind' in system:
		instance_numanode_membind = system['instance_numanode_membind']
		if isinstance(instance_numanode_membind, dict):
			instance_numanode_membind = instance_numanode_membind[cores_per_cpu]
		if isinstance(instance_numanode_membind, list):
			instance_numanode_membind = instance_numanode_membind[instance]

	# Set cpu_core to first core
	cpu_core = cores_per_cpu*instance_numanode_bind if isinstance(instance_numanode_bind, int) else 0
	if isinstance(system_ccc_config, dict) and 'instance_cpu_core_0' in system_ccc_config:
		if instance >= len(system_ccc_config['instance_cpu_core_0']):
			print('{} only defines {} instances for system {} ({} core) in {}'.format('instance_cpu_core_0', len(system_ccc_config['instance_cpu_core_0']), system_name, cores_per_cpu, config_filename))
			exit(1)
		
		cpu_core = system_ccc_config['instance_cpu_core_0'][instance]
	print('Fallback thread-masks start from CPU core {}.'.format(cpu_core))

	# Gather the list of threads for the instance
	if subsystem != '':
		assert 'subsystem_threads' in system, '\'subsystem_threads\' not defined for root-system {} in {}'.format(system_name, config_filename)
		assert subsystem in system['subsystem_threads'], '\'{}\' not defined in the subsystem_threads for root-system {} in {}'.format(subsystem, system_name, config_filename)
		threads = system['subsystem_threads'][subsystem]
	else:
		assert 'threads' in system, '\'threads\' not defined for system {} in {}'.format(system_name, config_filename)
		threads = system['threads']

	# Gather dictionaries of thread masks' lengths and user-specified thread masks
	thread_mask_length_dict = system['thread_mask_lengths'] if 'thread_mask_lengths' in system else {}
	thread_instance_mask_dict = None
	thread_instance_mask_dict_name = None
	if isinstance(system_ccc_config, dict):
		if subsystem == '' and 'instance_thread_masks' in system_ccc_config:
			thread_instance_mask_dict = system_ccc_config['instance_thread_masks']
			thread_instance_mask_dict_name = 'instance_thread_masks'
		elif 'instance_subsystem_thread_masks' in system_ccc_config and subsystem in system_ccc_config['instance_subsystem_thread_masks']:
			thread_instance_mask_dict = system_ccc_config['instance_subsystem_thread_masks'][subsystem]
			thread_instance_mask_dict_name = 'instance_subsystem_thread_masks'

	# Create the thread-list segment of the instance command
	hpguppi_threads_cmd_segment = []
	for thread in threads:
		if thread_instance_mask_dict is not None and (thread_instance_mask_dict is False or thread in thread_instance_mask_dict): # explicit core masks
			if thread_instance_mask_dict is False:
				if thread == threads[0]:
					print('{}[{}] defined as false: no thread-masking performed for system {} ({} core) in {}.'.format(
																					thread_instance_mask_dict_name, thread, system_full_name, cores_per_cpu, config_filename
																					)
																				)
				hpguppi_threads_cmd_segment.append(str(thread))
			else:
				if thread_instance_mask_dict[thread] is False:
					hpguppi_threads_cmd_segment.append(thread)
					continue
				
				assert isinstance(thread_instance_mask_dict[thread], list), '{}[{}] must define a mask for each instance as a list for system {} ({} core) in {}.'.format(
																					thread_instance_mask_dict_name, thread, system_full_name, cores_per_cpu, config_filename
																				)
				assert instance < len(thread_instance_mask_dict[thread]), '{}[{}] doesn\'t define a mask for instance {} for system {} ({} core) in {}.'.format(
																					thread_instance_mask_dict_name, thread, instance, system_full_name, cores_per_cpu, config_filename
																				)
				thread_mask = thread_instance_mask_dict[thread][instance]
				if isinstance(thread_mask, int):
					if thread not in thread_mask_length_dict:
						hpguppi_threads_cmd_segment.extend(['-c', str(thread_mask), str(thread)])
						cpu_core = thread_mask + 1
					else: # implies a mask of 'thread_mask_length' consequetive cores
						mask_val = 0
						for core_idx in range(thread_mask, thread_mask+thread_mask_length_dict[thread]):
							mask_val += 2**core_idx
						hpguppi_threads_cmd_segment.extend(['-m', str(mask_val), str(thread)])
						cpu_core = thread_mask + thread_mask_length_dict[thread]
				elif isinstance(thread_mask, list):
					mask_val = 0
					for core_idx in thread_mask:
						mask_val += 2**core_idx
						cpu_core = core_idx + 1
					hpguppi_threads_cmd_segment.extend(['-m', str(mask_val), str(thread)])
		else: # sequential core masks
			if thread_instance_mask_dict is not None:
				print('Fallback masking thread {} from core {}'.format(thread, cpu_core))
			if thread not in thread_mask_length_dict:
				hpguppi_threads_cmd_segment.extend(['-c', str(cpu_core), str(thread)])
				cpu_core += 1
			else:
				mask_val = 0
				for core_idx in range(cpu_core, cpu_core+thread_mask_length_dict[thread]):
					mask_val += 2**core_idx
				hpguppi_threads_cmd_segment.extend(['-m', str(mask_val), str(thread)])
				cpu_core += thread_mask_length_dict[thread]

	# Gather instance-agnostic instantiation variables
	prefix_exec = system.get('prefix_exec', default_prefix_exec)
	prefix_lib = system.get('prefix_lib', default_prefix_lib)
	logdir = system['logdir'] if 'logdir' in system else default_logdir
	# Gather optional instance-agnostic instantiation variables
	command_prefix = system.get('command_prefix', '')
	hpguppi_plugin = system.get('hpguppi_plugin', 'hpguppi_daq.so')

	# Gather required instance-sensitive instantiation variables
	assert 'instance_datadir' in system, '{} for system {} ({} core) in {}'.format('instance_datadir', system_name, cpu_core_count, config_filename)
	instance_datadir = system['instance_datadir']
	if isinstance(instance_datadir, dict):
		instance_datadir = instance_datadir[cores_per_cpu]
	instance_datadir = instance_datadir[instance]

	assert os.path.exists(instance_datadir), '{} datadir path does not exist for instance {} of system {} ({} core) in {}'.format(
		instance_datadir, instance, system_name, cores_per_cpu, config_filename)

	assert 'instance_bindhost' in system, '{} for system {} ({} core) in {}'.format('instance_bindhost', system_name, cpu_core_count, config_filename)
	instance_bindhost = system['instance_bindhost']
	if isinstance(instance_bindhost, dict):
		instance_bindhost = instance_bindhost[cores_per_cpu]
	instance_bindhost = instance_bindhost[instance]


	# Gather optional instance-sensitive instantiation variables
	instance_port_bind = None
	if 'instance_port_bind' in system:
		instance_port_bind = system['instance_port_bind'][instance]
	else:
		print('{} not found for system {} in {}, no per-instance BINDPORT option set'.format('instance_port_bind', system_name, config_filename))

	# Build options (key=value pairs for the status buffer)
	options = [
		'DATADIR={}'.format(instance_datadir),
		'BINDHOST={}'.format(instance_bindhost),
	]
	if instance_port_bind is not None:
		options.append('BINDPORT={}'.format(instance_port_bind))

	_keyword_variable_dict = {
		'BINDHOST': instance_bindhost,
		'INSTANCE': instance,
	}

	if 'options' in system:
		options += collect_yaml_value(system['options'], instance, cores_per_cpu, subsystem)

		# replace keywords
		for (option_idx, option) in enumerate(options):
			for var,keyword_val in _keyword_variable_dict.items():
					option = option.replace('${}'.format(var), str(keyword_val))
			options[option_idx] = option

	# Print empty line to conclude setup and assumption prints
	print()

	# Build hpguppi_daq command
	cmd = [
		command_prefix,
		'{}hashpipe'.format(prefix_exec),
		'-p',
		os.path.join(prefix_lib, hpguppi_plugin),
		'-I',
		str(instance),
	]
	
	for opt in options:
		cmd.extend(['-o', '{}'.format(opt)])
	cmd.extend(additional_arguments)
	cmd.extend(hpguppi_threads_cmd_segment)

	if isinstance(instance_numanode_bind, str):
		cmd = ['numactl', '--interleave={}'.format(instance_numanode_bind)] + cmd
	elif instance_numanode_bind is not False:
		cmd = ['numactl', '--cpunodebind={}'.format(instance_numanode_bind), '--membind={}'.format(instance_numanode_membind)] + cmd

	cmd = [seg for seg in cmd if seg != '']

	# Kill previous instances
	previous_instance_cmd_pattern = 'hashpipe -p .*\.so -I {}'.format(instance)
	if not dry_run:
		print(subprocess.run(['pkill', '-ef', previous_instance_cmd_pattern], capture_output=True).stdout.decode())

	# Handle logging true switch
	out_logpath = None
	err_logpath = None
	if logdir is not None:
		# Generate log filepaths
		hostname = socket.gethostname()
		out_logpath = os.path.join(logdir, '{}.{}.out'.format(hostname, instance))
		err_logpath = os.path.join(logdir, '{}.{}.err'.format(hostname, instance))

		# Trim existing logs
		for logpath in [out_logpath, err_logpath]:
			if os.path.exists(logpath):
				print('Trimming log: {}'.format(logpath))

				with open('/tmp/tmp.out', 'w') as tmpio:
					subprocess.run('tail -n 100000 {}'.format(logpath).split(' '), stdout=tmpio)
				
				subprocess.run('mv /tmp/tmp.out {}'.format(logpath).split(' '))
				with open(logpath, 'a') as logio:
					logio.write('\n{}\nStartup {}\n{}\n{}\n'.format('-'*20, datetime.datetime.now(), ' '.join(cmd), 'v'*20))
		print()
	
	out_logio = None if out_logpath is None else open(out_logpath, 'a')
	err_logio = None if err_logpath is None else open(err_logpath, 'a')

	# Setup environment
	environment_keys = environment_keys
	system_environment_keys_start_index = len(environment_keys)
	if 'hashpipe_keyfile' in system:
		environment_keys.append('HASHPIPE_KEYFILE={}'.format(system['hashpipe_keyfile']))
	if 'environment' in system:
		environment_keys += collect_yaml_value(system['environment'], instance, cores_per_cpu, subsystem)

	hashpipe_env = get_extended_environment(environment_keys, _keyword_variable_dict)

	if 'setup_commands' in system:
		for setup_command in collect_yaml_value(system['setup_commands'], instance, cores_per_cpu, subsystem):
			for var,val in _keyword_variable_dict.items():
				setup_command = setup_command.replace('${}'.format(var), str(val))
			print('#', setup_command)
			out_logio.write('\n# {}\n'.format(setup_command))
			if not dry_run:
				out_logio.write(subprocess.run(setup_command.split(' '), env=hashpipe_env, capture_output=True).stdout.decode())
		print()
		out_logio.write('%'*20+'\n')

	print(' '.join(environment_keys[system_environment_keys_start_index:] + cmd))

	if not dry_run:
		subprocess.Popen(cmd, env=hashpipe_env, stdout=out_logio, stderr=err_logio)
	else:
		print(hashpipe_env)
		print('^^^ Dry run ^^^')
		err_logio.write('Dry run')
		out_logio.write('Dry run')
		err_logio.close()
		out_logio.close()

	if 'post_commands' in system:
		post_commands = system['post_commands']
		post_command_index = 0
		while post_command_index < len(post_commands):
			post_command = post_commands[post_command_index]
			post_command_index += 1

			if isinstance(post_command, dict): # cores_per_cpu dict, insert commands next and continue
				assert cores_per_cpu in post_command, 'Missing an entry for {} cores in the {} Dict of for system {} in {}\n'.format(cores_per_cpu, 'post_commands', system_name, config_filename, post_command)
				if isinstance(post_command[cores_per_cpu], list):
					post_commands[post_command_index:post_command_index] = post_command[cores_per_cpu]
				else:
					post_commands.insert(post_command_index, post_command[cores_per_cpu])
				continue
			if isinstance(post_command, list): # instance-indexed list of commands, select appropriately
				post_command = post_command[instance]
			
			for var,val in _keyword_variable_dict.items():
				post_command = post_command.replace('${}'.format(var), str(val))
			print('#', post_command)
			if not dry_run:
				print(subprocess.run(post_command.split(' '), env=hashpipe_env, capture_output=True).stdout.decode())
		print()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Starts an instance of Hpguppi_daq')
	parser.add_argument('system', type=str,
											help='The name of the system the instance will be configured for')
	parser.add_argument('instance', type=str,
											help='The instance enumeration and ID (\',\' delimited)')
	parser.add_argument('additional_arguments', type=str, nargs='*',
											help='Arguments added to the instance command')
	parser.add_argument('-o', '--options', type=str, nargs='*', default=[],
											help='Additional key=val pairs for the Hpguppi instance\'s status buffer')
	parser.add_argument('-e', '--environment-keys', type=str, nargs='*', default=[],
											help='Additional key=val pairs for the environment of the Hpguppi instance')
	parser.add_argument('-d', '--dry-run', action='store_true',
											help='Do not start the instance')
	parser.add_argument('-D', '--delete', action='store_true',
											help='Delete the shared memory of the instance')
	parser.add_argument('--configfile', type=str, default='config_hpguppi.yml',
											help='The file containing systems\' configurations.')
	args = parser.parse_args()

	
	# Load configuration file 
	with open(args.configfile, 'r') as fio:
		yaml_config = yaml.load(fio, Loader=yaml.SafeLoader)
		for inst_str in args.instance.split(','):
			if args.delete:
				delete(
					args.system,
					int(inst_str),
					yaml_config,
					environment_keys = args.environment_keys,
					dry_run = args.dry_run,
					config_filename = args.configfile
				)
			else:
				run(
					args.system,
					int(inst_str),
					yaml_config,
					additional_arguments = args.additional_arguments,
					options = args.options,
					environment_keys = args.environment_keys,
					dry_run = args.dry_run,
					config_filename = args.configfile
				)

	exit(0)