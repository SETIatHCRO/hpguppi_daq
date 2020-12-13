import subprocess
import argparse
import time
import os
import glob
import re

def get_hashpipe_key_value(key):
	return subprocess.run(['hashpipe_check_status', '--query='+key], capture_output=True).stdout

def get_hashpipe_pulse():
	return get_hashpipe_key_value('DAQPULSE')

def block_until_pulse_change(maxstale=20, silent=False):
	stalecount = 0
	stalepulse = get_hashpipe_pulse()
	while True:
		try:
			stalepulse.decode()
			break
		except:
			stalecount += 1
		if not silent:
			print('\rDAQPULSE not a UTF-8 string'+'.'*(stalecount%5), end=' '*5)
		time.sleep(2)
		if stalecount > 10:
			if not silent:
				print('\rExiting after excessive waiting for DAQPULSE to be set. (20 seconds)')
			return False
		stalepulse = get_hashpipe_pulse()

	stalecount = 0

	while get_hashpipe_pulse() == stalepulse:
		if not silent:
			print('\rwaiting for DAQPULSE to change'+'.'*(stalecount%5), end=' '*5)
		time.sleep(1)
		if (stalecount > maxstale):
			if not silent:
				print('\rExiting after excessive waiting for DAQPULSE to change. (%d seconds)' %(maxstale))
			break
		stalecount += 1
	if not silent:
		print('')
	if (stalecount > maxstale):
		return False
	return True

def start_hashpipe(bindhost):
	"""Does not directly kill existing hashpipes"""
	print('\nStarting hashpipe (', bindhost, ')', sep='')
	subprocess.run(['sudo', '../src/init_atasnap_hashpipe.sh', '-o', 'BINDHOST='+bindhost]) 

def kill_hashpipes():
	print('\nKilling existing hashpipes')
	subprocess.run(['sudo', 'pkill', '-e', '-f', '.*/hashpipe\s.*'])

def start_redis_gateway():
	subprocess.run(['../src/init_redis_gateway.sh']) # assume it will kill existing gateways

def kill_hashpipe_relevant():
	print('\nKilling existing hashpipe-related processes')
	subprocess.run(['sudo', 'pkill', '-e', '-f', '.*hashpipe.*'])

def clear_shared_memory():
	print('\nDeleting shared memory')
	lines = subprocess.run(['ipcs'], capture_output=True).stdout.decode().split('\n')
	mode = 0
	for line in lines:
		if mode == 1 or mode == 3:
			m = re.match(r'[^\s]*\s+([^\s]+).*', line)
			if m:
				flag = '-m' if mode == 1 else '-s'
				subprocess.run(['sudo', 'ipcrm', flag, str(m.group(1))])
			else:
				mode += 1
		elif mode == 0 and re.match(r'.*shmid.*', line):
			mode = 1
		elif mode == 2 and re.match(r'.*semid.*', line):
			mode = 3
		elif mode == 4:
			break

def get_hashpipe_capture_dir():
	rawfiledir = ''
	for key in ['DATADIR', 'PROJID', 'BACKEND', 'BANKNAM']:
		part = get_hashpipe_key_value(key)
		try:
			rawfiledir = os.path.join(rawfiledir, part.decode().strip())
		except:
			print('Status Key', key, 'is not a valid string. Exiting.')
			rawfiledir = False
	return rawfiledir

def get_latest_raw_stem_in_dir(rawfiledir):
	files = glob.glob(rawfiledir+'/*.raw')
	files.sort(key=os.path.getmtime, reverse=True)

	if len(files) == 0:
		print('Could not find any *.raw files in "', rawfiledir, '".', sep="")
		exit(1)

	latest = files[0]
	stemmatch = re.match(r'(.*)\.\d{4}\.raw', latest, re.M|re.I)
	if not stemmatch:
		print('Could not match a stem pattern against "', latest, '".', sep="")
		return False
	return stemmatch.group(1)
