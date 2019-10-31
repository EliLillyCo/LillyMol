#include <memory>
#include <limits>

#include "accumulator.h"
#include "misc.h"

#include "substructure.h"
#include "spatially_common_matched_atoms.h"

static double
proximity_to_most(const Molecule & m,
                   const Set_of_Atoms & e,
                   Accumulator<double> & x,
                   Accumulator<double> & y,
                   Accumulator<double> & z)
{
  const int n = e.number_elements();

  double rc = 0.0;

  const float xmean = x.average();
  const float ymean = y.average();
  const float zmean = z.average();

  for (int i = 0; i < n; ++i)
  {
    const Atom * a = m.atomi(e[i]);

    rc += (a->x() - xmean) * (a->x() - xmean) + (a->y() - ymean) * (a->y() - ymean) + (a->z() - zmean) * (a->z() - zmean);
  }

  return rc;
}

static int
identify_embedding_closest_to_most(const Molecule & m,
                   Substructure_Results & sresults,
                   Accumulator<double> & x,
                   Accumulator<double> & y,
                   Accumulator<double> & z)
{
  const int n = sresults.number_embeddings();

  double min_dist = std::numeric_limits<double>::max();
  int id_with_min_dist = -1;

  for (int i = 0; i < n; ++i)
  {
    const Set_of_Atoms * e = sresults.embedding(i);

    const double d = proximity_to_most(m, *e, x, y, z);

    if (d < min_dist)
    {
      min_dist = d;
      id_with_min_dist = i;
    }
  }

//std::cerr << "Start with " << n << " embeddings, min dist " << min_dist << " for embedding " << id_with_min_dist << std::endl;

  for (int i = 0; i < id_with_min_dist; ++i)
  {
    sresults.remove_embedding(0);
  }
  while (sresults.number_embeddings() > 1)
  {
    sresults.remove_embedding(1);
  }

  return 1;
}

static void
gather_coordinates(const Molecule & m,
                   const Set_of_Atoms & e,
                   Accumulator<double> & x,
                   Accumulator<double> & y,
                   Accumulator<double> & z)
{
  const int n = e.number_elements();

  for (int i = 0; i < n; ++i)
  {
    const Atom * a = m.atomi(e[i]);

    x.extra(a->x());
    y.extra(a->y());
    z.extra(a->z());
  }

  return;
}

int
spatially_common_matched_atoms(resizable_array_p<Molecule> & molecules,
                        Substructure_Results * sresults)
{
  const int n = molecules.number_elements();

  int * multiple_embeddings = new_int(n); std::unique_ptr<int[]> free_multiple_embeddings(multiple_embeddings);

  int max_embeddings = 0;

  for (int i = 0; i < n; ++i)
  {
    multiple_embeddings[i] = sresults[i].number_embeddings();
    if (multiple_embeddings[i] > max_embeddings)
      max_embeddings = multiple_embeddings[i];
  }

  if (1 == max_embeddings)
    return 0;

  Accumulator<double> x, y, z;

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < multiple_embeddings[i]; ++j)
    {
      const Set_of_Atoms * e = sresults[i].embedding(j);

      gather_coordinates(*molecules[i], *e, x, y, z);
    }
  }

  int rc = 0;

  for (int i = 0; i < n; ++i)
  {
    if (1 == multiple_embeddings[i])
      continue;

    identify_embedding_closest_to_most(*molecules[i], sresults[i], x, y, z);

    rc++;
  }

  return rc;
}
