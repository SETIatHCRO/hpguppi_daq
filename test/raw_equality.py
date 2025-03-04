from blimpy import GuppiRaw
import numpy as np
import argparse

def read_datablock(guppi: GuppiRaw):
  header, dx, dy = guppi.read_next_data_block_int8()
  if dx is None:
    return None, None 
  return header, np.append(dx, dy, axis=2)

def read_data(guppi: GuppiRaw):
  _, data = read_datablock(guppi)
  return data

def next_data_isequal(ref_gr: GuppiRaw, oth_gr: GuppiRaw):
  """
  Compares the next section of data between 2 RAW files, accounting
  for different BLOCSIZEs. For the last comparison `None` is returned
  if both files are completed, otherwise `False` is returned.
  """

  ref_header, ref_data = read_datablock(ref_gr)
  oth_header, oth_data = read_datablock(oth_gr)
  if ref_data is None and oth_data is None:
    return None
  elif ref_data is None or oth_data is None:
    return False

  ref_blksize = ref_header['BLOCSIZE']
  oth_blksize = oth_header['BLOCSIZE']

  if ref_blksize > oth_blksize:
    ratio, modulo = divmod(ref_blksize, oth_blksize)
    assert modulo == 0, f'BLOCSIZE values are not whole multiples: `{ref_blksize}`/`{oth_blksize}`'
    ntime = oth_data.shape[1]
    first_equal = np.equal(ref_data[:, 0:ntime, :], oth_data).all()
    return first_equal and all(np.equal(ref_data[:, i*ntime:(i+1)*ntime, :], read_data(oth_gr)).all() for i in range(1, ratio))
  elif ref_blksize < oth_blksize:
    ratio, modulo = divmod(oth_blksize, ref_blksize)
    assert modulo == 0, f'BLOCSIZE values are not whole multiples: `{oth_blksize}`/`{ref_blksize}`'
    ntime = ref_data.shape[1]
    first_equal = np.equal(oth_data[:, 0:ntime, :], ref_data).all()
    return first_equal and all(np.equal(oth_data[:, i*ntime:(i+1)*ntime, :], read_data(ref_gr)).all() for i in range(1, ratio))
  else:
    return np.equal(ref_data == oth_data).all()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Validates the data between 2 GUPPI RAW files.')
  parser.add_argument('reference_filepath', type=str,
                      help='The filepath of the reference GUPPI RAW file.')
  parser.add_argument('test_filepath', type=str,
                      help='The filepath of the GUPPI RAW file to be validated.')
  args = parser.parse_args()

  print(f'Opening `{args.reference_filepath}`')
  ref_gr = GuppiRaw(args.reference_filepath)
  print(f'Opening `{args.test_filepath}`')
  test_gr = GuppiRaw(args.test_filepath)

  block_num=0
  while True:
    equality = next_data_isequal(ref_gr, test_gr)

    if equality is None:
      break
    elif equality is False:
      print(f'Block #{block_num} is not equivalent.')
      exit(1)
    # print(f'Block #{block_num} is equivalent.')
    block_num += 1

  print(f'{block_num} blocks verified as equal.')
  exit(0)