#include <stdlib.h>

#include "Foundational/iwmisc/misc.h"

#include "gfp.h"

/*
  We are going to copy a gfp fingerprint to contiguous storage.
  We need 5 words to hold
    molecular properties integer present or absent
    molecular properties continuous present or absent
    number fixed size fingerprints
    number sparse fingerprints
*/

int
IW_General_Fingerprint::bytes_needed_for_contiguous_storage() const
{
  int rc = 4 * sizeof(int);

  if (_molecular_properties_integer.active())
    rc += sizeof(Molecular_Properties<int>) + _molecular_properties_integer.nproperties() * sizeof(int);

  if (_molecular_properties_continuous.active())
    rc += sizeof(Molecular_Properties<float>) + _molecular_properties_continuous.nproperties() * sizeof(float);

  int n = number_fingerprints();

  for (int i = 0; i < n; i++)
  {
    rc += sizeof(IWDYFP) + _fingerprint[i].nbits() / IW_BITS_PER_BYTE;
  }

  n = number_sparse_fingerprints();

  for (int i = 0; i < n; i++)
  {
    int nb = _sparse_fingerprint[i].nbits();

    rc += sizeof(Sparse_Fingerprint);

    rc += nb * 2 * sizeof(int);    // _bit and _count are each arrays
  }

  return rc;
}

void *
IW_General_Fingerprint::copy_to_contiguous_storage (void * p) const
{
  int * iptr = reinterpret_cast<int *>(p);

  int npri = _molecular_properties_integer.nproperties();

  *iptr = npri;
  iptr++;

  int nprf = _molecular_properties_continuous.nproperties();
  *iptr = nprf;
  iptr++;

  *iptr = number_fingerprints();
  iptr++;

  *iptr = number_sparse_fingerprints();
  iptr++;

  if (npri)
    iptr = (int *) _molecular_properties_integer.copy_to_contiguous_storage(iptr);

  if (nprf)
    iptr = (int *) _molecular_properties_continuous.copy_to_contiguous_storage(iptr);

  int n = number_fingerprints();

  for (int i = 0; i < n; i++)
  {
    iptr = (int *) _fingerprint[i].copy_to_contiguous_storage(iptr);
  }

  n = number_sparse_fingerprints();
  for (int i = 0; i < n; i++)
  {
    iptr = (int *) _sparse_fingerprint[i].copy_to_contiguous_storage(iptr);
  }

  return iptr;
}

int
IW_General_Fingerprint::bytes_needed_for_contiguous_storage_gpu() const
{
  int rc = 4 * sizeof(int);

  if (_molecular_properties_integer.active())
    rc += sizeof(int) + _molecular_properties_integer.nproperties() * sizeof(int);

  if (_molecular_properties_continuous.active())
    rc += sizeof(int) + _molecular_properties_continuous.nproperties() * sizeof(float);

  int n = number_fingerprints();

  for (int i = 0; i < n; i++)
  {
    rc += sizeof(int) + _fingerprint[i].nbits() / IW_BITS_PER_BYTE + sizeof(int);
  }

  n = number_sparse_fingerprints();

  for (int i = 0; i < n; i++)
  {
    int nb = _sparse_fingerprint[i].nbits();

    rc += 3*sizeof(int);

    rc += nb * 2 * sizeof(int);    // _bit and _count are each arrays
  }

  return rc;
}

void *
IW_General_Fingerprint::copy_to_contiguous_storage_gpu (void * p) const
{
  int * iptr = reinterpret_cast<int *>(p);

  int npri = _molecular_properties_integer.nproperties();

  *iptr = npri;
  iptr++;

  int nprf = _molecular_properties_continuous.nproperties();
  *iptr = nprf;
  iptr++;

  *iptr = number_fingerprints();
  iptr++;

  *iptr = number_sparse_fingerprints();
  iptr++;

  if (npri)
    iptr = (int *) _molecular_properties_integer.copy_to_contiguous_storage_gpu(iptr);

  if (nprf)
    iptr = (int *) _molecular_properties_continuous.copy_to_contiguous_storage_gpu(iptr);

  int n = number_fingerprints();

  for (int i = 0; i < n; i++)
  {
    iptr = (int *) _fingerprint[i].copy_to_contiguous_storage_gpu(iptr);
  }

  n = number_sparse_fingerprints();
  for (int i = 0; i < n; i++)
  {
    iptr = (int *) _sparse_fingerprint[i].copy_to_contiguous_storage_gpu(iptr);
  }

  return iptr;
}

const void *
IW_General_Fingerprint::build_from_contiguous_storage (const void * p,
                                        int allocate_arrays)
{
  const int * iptr = reinterpret_cast<const int *>(p);

  int npri = *iptr;
  iptr++;
  int nprf = *iptr;
  iptr++;
  int nfixed = *iptr;
  iptr++;
  int nsparse = *iptr;
  *iptr++;

//cerr << "Copying fingerprint with " << npri << " integer properties, " << nfixed << " fixed and " << nsparse << " sparse fingerprints\n";

  if (npri)
    iptr = (int *) _molecular_properties_integer.build_from_contiguous_storage(iptr, allocate_arrays);

  if (nprf)
    iptr = (int *) _molecular_properties_continuous.build_from_contiguous_storage(iptr, allocate_arrays);

  for (int i = 0; i < nfixed; i++)
  {
    iptr = (int *) _fingerprint[i].build_from_contiguous_storage(iptr, allocate_arrays);
  }

  for (int i = 0; i < nsparse; i++)
  {
    iptr = (int *) _sparse_fingerprint[i].build_from_contiguous_storage(iptr, allocate_arrays);
  }

  return iptr;
}
