/* hpguppi_databuf.c
 *
 * Routines for creating and accessing main data transfer
 * buffer in shared memory.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <errno.h>
#include <time.h>

#include <hashpipe.h>

#include "hpguppi_databuf.h"
//#include "fitshead.h"

hashpipe_databuf_t *hpguppi_input_databuf_create(int instance_id, int databuf_id)
{
    int i;

    /* Calc databuf sizes */
    size_t header_size = sizeof(hashpipe_databuf_t)
                       + sizeof(hashpipe_databuf_alignment);
    size_t block_size  = sizeof(hpguppi_input_block_t);
    int    n_block = N_INPUT_BLOCKS;

    hpguppi_input_databuf_t * d = (hpguppi_input_databuf_t *)
        hashpipe_databuf_create(
            instance_id, databuf_id, header_size, block_size, n_block);

    if(!d) {
      return NULL;
    }

    /* Zero out blocks */
    for(i=0; i<n_block; i++) {
      memset(&(d->block[i]), 0, sizeof(hpguppi_input_block_t));
    }

    /* Init headers of each databuf block */
    char end_key[81];
    memset(end_key, ' ', 80);
    strncpy(end_key, "END", 3);
    end_key[80]='\0';
    for (i=0; i<n_block; i++) { 
        memcpy(hpguppi_databuf_header(d,i), end_key, 80); 
    }

    return (hashpipe_databuf_t *)d;
}

hashpipe_databuf_t *hpguppi_databuf_attach_retry(int instance_id, int databuf_id) {
  // Attach to databuf as a low-level hashpipe databuf. Cannot create
  // the upstream databuf if it does not yet exist. Wait at most 1 second
  // for it to be created by a different thread. 
  struct timespec ts = {0, 1000}; // One microsecond
  int max_tries = 1000000; // One million microseconds
  hashpipe_databuf_t *db = NULL;
  for(int i = 0; i < max_tries; i++) {
      db = hashpipe_databuf_attach(instance_id, databuf_id);
      if(db) break;
      nanosleep(&ts, NULL);
  }

  return db;
}

#if 0 // OLD STUFF

int hpguppi_databuf_detach(struct guppi_databuf *d) {
    int rv = shmdt(d);
    if (rv!=0) {
        guppi_error("guppi_status_detach", "shmdt error");
        return(GUPPI_ERR_SYS);
    }
    return(GUPPI_OK);
}

void guppi_databuf_clear(struct guppi_databuf *d) {

    /* Zero out semaphores */
    union semun arg;
    arg.array = (unsigned short *)malloc(sizeof(unsigned short)*d->n_block);
    memset(arg.array, 0, sizeof(unsigned short)*d->n_block);
    semctl(d->semid, 0, SETALL, arg);
    free(arg.array);

    /* Clear all headers */
    int i;
    for (i=0; i<d->n_block; i++) {
        guppi_fitsbuf_clear(guppi_databuf_header(d, i));
    }

}

void guppi_fitsbuf_clear(char *buf) {
    char *end, *ptr;
    end = ksearch(buf, "END");
    if (end!=NULL) {
        for (ptr=buf; ptr<=end; ptr+=80) memset(ptr, ' ', 80);
    }
    memset(buf, ' ' , 80);
    strncpy(buf, "END", 3);
}

char *guppi_databuf_header(struct guppi_databuf *d, int block_id) {
    return((char *)d + d->struct_size + block_id*d->header_size);
}

char *guppi_databuf_data(struct guppi_databuf *d, int block_id) {
    return((char *)d + d->struct_size + d->n_block*d->header_size
            + block_id*d->block_size);
}

struct guppi_databuf *guppi_databuf_attach(int databuf_id) {

    /* Get shmid */
    int shmid;
    shmid = shmget(GUPPI_DATABUF_KEY + databuf_id - 1, 0, 0666);
    if (shmid==-1) {
        // Doesn't exist, exit quietly otherwise complain
        if (errno!=ENOENT)
            guppi_error("guppi_databuf_attach", "shmget error");
        else
            guppi_error("guppi_databuf_attach", "does not exist");
        return(NULL);
    }

    /* Attach */
    struct guppi_databuf *d;
    d = shmat(shmid, NULL, 0);
    if (d==(void *)-1) {
        guppi_error("guppi_databuf_attach", "shmat error");
        return(NULL);
    }

    return(d);

}

int guppi_databuf_block_status(struct guppi_databuf *d, int block_id) {
    return(semctl(d->semid, block_id, GETVAL));
}

int guppi_databuf_total_status(struct guppi_databuf *d) {

    /* Get all values at once */
    union semun arg;
    arg.array = (unsigned short *)malloc(sizeof(unsigned short)*d->n_block);
    memset(arg.array, 0, sizeof(unsigned short)*d->n_block);
    semctl(d->semid, 0, GETALL, arg);
    int i,tot=0;
    for (i=0; i<d->n_block; i++) tot+=arg.array[i];
    free(arg.array);
    return(tot);

}

int guppi_databuf_bitmask_status(struct guppi_databuf *d) {

    /* Get all values at once */
    union semun arg;
    arg.array = (unsigned short *)malloc(sizeof(unsigned short)*d->n_block);
    memset(arg.array, 0, sizeof(unsigned short)*d->n_block);
    semctl(d->semid, 0, GETALL, arg);
    int i,bitmask=0;
    for (i=0; i<d->n_block && i<8*sizeof(bitmask); i++) {
      if(arg.array[i]) {
        bitmask |= (1<<i);
      }
    }
    free(arg.array);
    return bitmask;
}

void guppi_databuf_str_status(struct guppi_databuf *d, char * str, size_t len) {
    /* Get all values at once */
    union semun arg;
    arg.array = (unsigned short *)malloc(sizeof(unsigned short)*d->n_block);
    memset(arg.array, 0, sizeof(unsigned short)*d->n_block);
    semctl(d->semid, 0, GETALL, arg);
    int i;
    for (i=0; i<d->n_block && i<len-1; i++) {
        str[i] = arg.array[i] ? 'X' : '.';
    }
    // Nul terminate
    str[i] = '\0';

    free(arg.array);
}

int guppi_databuf_wait_free(struct guppi_databuf *d, int block_id) {
    int rv;
    struct sembuf op;
    op.sem_num = block_id;
    op.sem_op = 0;
    op.sem_flg = 0;
    struct timespec timeout;
    timeout.tv_sec = 0;
    timeout.tv_nsec = 250000000;
    rv = semtimedop(d->semid, &op, 1, &timeout);
    if (rv==-1) { 
        if (errno==EAGAIN) return(GUPPI_TIMEOUT);
        if (errno==EINTR) return(GUPPI_ERR_SYS);
        guppi_error("guppi_databuf_wait_free", "semop error");
        perror("semop");
        return(GUPPI_ERR_SYS);
    }
    return(0);
}

int guppi_databuf_wait_filled(struct guppi_databuf *d, int block_id) {
    /* This needs to wait for the semval of the given block
     * to become > 0, but NOT immediately decrement it to 0.
     * Probably do this by giving an array of semops, since
     * (afaik) the whole array happens atomically:
     * step 1: wait for val=1 then decrement (semop=-1)
     * step 2: increment by 1 (semop=1)
     */
    int rv;
    struct sembuf op[2];
    op[0].sem_num = op[1].sem_num = block_id;
    op[0].sem_flg = op[1].sem_flg = 0;
    op[0].sem_op = -1;
    op[1].sem_op = 1;
    struct timespec timeout;
    timeout.tv_sec = 0;
    timeout.tv_nsec = 250000000;
    rv = semtimedop(d->semid, op, 2, &timeout);
    if (rv==-1) { 
        if (errno==EAGAIN) return(GUPPI_TIMEOUT);
        // Don't complain on a signal interruption
        if (errno==EINTR) return(GUPPI_ERR_SYS);
        guppi_error("guppi_databuf_wait_filled", "semop error");
        perror("semop");
        return(GUPPI_ERR_SYS);
    }
    return(0);
}

int guppi_databuf_set_free(struct guppi_databuf *d, int block_id) {
    /* This function should always succeed regardless of the current
     * state of the specified databuf.  So we use semctl (not semop) to set
     * the value to zero.
     */
    int rv;
    union semun arg;
    arg.val = 0;
    rv = semctl(d->semid, block_id, SETVAL, arg);
    if (rv==-1) { 
        guppi_error("guppi_databuf_set_free", "semctl error");
        return(GUPPI_ERR_SYS);
    }
    return(0);
}

int guppi_databuf_set_filled(struct guppi_databuf *d, int block_id) {
    /* This function should always succeed regardless of the current
     * state of the specified databuf.  So we use semctl (not semop) to set
     * the value to one.
     */
    int rv;
    union semun arg;
    arg.val = 1;
    rv = semctl(d->semid, block_id, SETVAL, arg);
    if (rv==-1) { 
        guppi_error("guppi_databuf_set_filled", "semctl error");
        return(GUPPI_ERR_SYS);
    }
    return(0);
}
#endif // OLD STUFF
