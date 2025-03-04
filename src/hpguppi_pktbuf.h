// hpguppi_ibverbs_pkt_thread.h
//
// Header file for externally visible data structures and functions in
// hpguppi_ibverbs_pkt_thread.c.

#ifndef _HPGUPPI_PKTBUF_H_
#define _HPGUPPI_PKTBUF_H_

#include <stdint.h>
#include <string.h>
#include <stddef.h>
#include <time.h>

#include "hashpipe.h"
#include "hpguppi_databuf.h"

// Alignment size to use.  Currently set to 64 (== 512/8) for compatibility
// with AVX512 instructions.
#define PKT_ALIGNMENT_SIZE (64)

// Maximum number of chunks supported
#define MAX_CHUNKS (8)

// Structure that holds info about a "chunk".  A chunk is part of a packet that
// is stored at a PKT_ALIGNMENT_SIZE aligned address.  The chunk_size is the
// number of bytes from the packet that are stored in the chunk.  The
// chunk_aligned_size is chunk_size rounded up to the next multple of
// PKT_ALIGNMENT_SIZE.  The chunk_offset is the offset of the chunk from the
// packet's first chunk.  The first chunk will have a chunk_offset of 0.
struct hpguppi_pktbuf_chunk {
  size_t chunk_size;
  size_t chunk_aligned_size;
  off_t chunk_offset;
};

// Structure that holds info about packet/slot/block sizing.  A block is
// divided into "slots".  Each slot holds one packet, possibly with internal
// and/or trailing padding added to align various sections of the packet to
// PKT_ALIGNMENT_SIZE.  These sections are called chunks.  num_chunks specifies
// the number of chunks that are being used.  The pkt_size is the sum of the
// (unaligned) sizes of all chunks.  The slot_size is the size of a slot and
// equals the sum of the chunk_aligned_sizes.  slots_per_block is the number of
// slots in a data block.  Note that slot_size * slots_per_block may be less
// than the size of data block by up PKT_ALIGNMENT_SIZE-1 bytes.
struct hpguppi_pktbuf_info {
  uint32_t num_chunks;
  size_t pkt_size;
  size_t slot_size;
  size_t slots_per_block;
  struct hpguppi_pktbuf_chunk chunks[MAX_CHUNKS];
};

// Function to get a pointer to a databuf's pktbuf_info structure.  Assumes
// that the pktbuf_info structure is tucked into the "padding" bytes of the
// hpguppi_intput_databuf.
// TODO Check db->header.data_type?
static inline
struct hpguppi_pktbuf_info *
hpguppi_pktbuf_info_ptr(hpguppi_input_databuf_t *db)
{
  return (struct hpguppi_pktbuf_info *)(db->padding);
}

// Function to get the offset within a slot to an (unaligned) offset within a
// packet.  This accounts for the padding between chunks.  For example, if the
// chuck sizes are 14,20,1500 (e.g. MAC,IP,PAYLOAD) and the chunks are aligned
// on 64 byte boundaries, then (unaligned) packet offset 34 (e.g. start of
// PAYLOAD) would have (aligned) slot offset 64.
off_t
hpguppi_pktbuf_slot_offset(hpguppi_input_databuf_t *db, off_t pkt_offset);

// Function that threads can call to wait for hpguppi_ibvpkt_thread to finalize
// ibverbs setup.  After this funtion returns, the underlying
// hashpipe_ibv_context structure used by hpguppi_ibvpkt_thread will be fully
// initialized and flows can be created/destroyed by calling
// hpguppi_ibvpkt_flow().
void hpguppi_ibvpkt_wait_running(hashpipe_status_t * st);

// Function to get a pointer to slot "slot_id" in block "block_id" of databuf
// "db".
static inline
uint8_t *
hpguppi_pktbuf_block_slot_ptr(hpguppi_input_databuf_t *db,
    uint64_t block_id, uint32_t slot_id)
{
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(db);
  block_id %= db->header.n_block;
  return (uint8_t *)db->block[block_id].data + slot_id * pktbuf_info->slot_size;
}

#endif // _HPGUPPI_PKTBUF_H_
