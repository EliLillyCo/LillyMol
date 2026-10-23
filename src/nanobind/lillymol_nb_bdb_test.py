import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(__file__))


def _runfile_path(path):
    runfiles_dir = os.environ.get("RUNFILES_DIR")
    if runfiles_dir:
        candidate = os.path.join(runfiles_dir, "_main", path)
        if os.path.exists(candidate):
            return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", path))
    if os.path.exists(candidate):
        return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bazel-bin", path))
    if os.path.exists(candidate):
        return candidate

    return None

import lillymol
import lillymol_bdb


class TestNanobindBdb(unittest.TestCase):

    def test_ring_database(self):
        tool = _runfile_path("Molecule_Tools_Bdb/smi2rings_bdb")
        self.assertIsNotNone(tool)
        with tempfile.TemporaryDirectory() as tmpdir:
            def build(name, smiles):
                source = os.path.join(tmpdir, name + ".smi")
                database = os.path.join(tmpdir, name + ".bdb")
                with open(source, "w", encoding="utf-8") as writer:
                    for i, smi in enumerate(smiles):
                        writer.write(f"{smi} mol{i}\n")
                subprocess.run([tool, "-d", database, "-d", "STORE", "-j", "ring",
                                "-n", source], check=True, capture_output=True)
                return database

            first = build("first", ["c1ccccc1"] * 2 + ["C1CCCCC1"] * 3)
            second = build("second", ["c1ccccc1"] * 4)
            benzene = lillymol.MolFromSmiles("c1ccccc1 benzene")
            mixed = lillymol.MolFromSmiles("c1ccccc1.C1CCCCC1 mixed")
            unseen = lillymol.MolFromSmiles("C1CCCC1 unseen")
            acyclic = lillymol.MolFromSmiles("CCO ethanol")
            original = (mixed.unique_smiles(), mixed.name())
            with lillymol_bdb.RingDatabase(first, bare_ring_fallback=True) as db:
                self.assertEqual(db.rarest_ring_count(benzene), 2)
                self.assertEqual(db.rarest_ring_count(mixed), 2)
                self.assertEqual(db.rarest_ring_count(unseen), 0)
                self.assertIsNone(db.rarest_ring_count(acyclic))
                self.assertEqual((mixed.unique_smiles(), mixed.name()), original)
                self.assertEqual(db.rarest_ring_count(
                    lillymol.MolFromSmiles("c1ccccc1.C1CCCC1")), 0)
            db.close()
            with self.assertRaises(RuntimeError):
                db.rarest_ring_count(benzene)
            with lillymol_bdb.RingDatabase([first, second], bare_ring_fallback=True) as db:
                self.assertEqual(db.rarest_ring_count(benzene), 4)
                self.assertEqual(db.rarest_ring_count(mixed), 3)
            with lillymol_bdb.RingDatabase(first, max_ring_size=5) as db:
                self.assertIsNone(db.rarest_ring_count(benzene))
            toluene = lillymol.MolFromSmiles("Cc1ccccc1 toluene")
            with lillymol_bdb.RingDatabase(first) as db:
                self.assertEqual(db.rarest_ring_count(toluene), 0)
            with lillymol_bdb.RingDatabase(first, bare_ring_fallback=True) as db:
                self.assertEqual(db.rarest_ring_count(toluene), 2)
            corrupt_input = os.path.join(tmpdir, "corrupt.txt")
            corrupt_db = os.path.join(tmpdir, "corrupt.bdb")
            with open(corrupt_input, "w", encoding="utf-8") as writer:
                writer.write(f"{benzene.unique_smiles()} not_a_ring_proto\n")
            subprocess.run([_runfile_path("BerkeleyDB/iwbdb_load"), "-d", corrupt_db,
                            corrupt_input], check=True, capture_output=True)
            with lillymol_bdb.RingDatabase(corrupt_db) as db:
                with self.assertRaises(RuntimeError):
                    db.rarest_ring_count(benzene)
            labelled = lillymol.MolFromSmiles("[13CH2]1CCCCC1 labelled")
            before = labelled.unique_smiles()
            with lillymol_bdb.RingDatabase(
                    first, remove_isotopes=True, attachment_labels=False,
                    bare_ring_fallback=True, standardize=True) as db:
                self.assertEqual(db.rarest_ring_count(labelled), 3)
                self.assertEqual(labelled.unique_smiles(), before)
            with lillymol_bdb.RingDatabase(first, standardize=["all"]) as db:
                self.assertIsNone(db.rarest_ring_count(acyclic))
            with self.assertRaises(ValueError):
                lillymol_bdb.RingDatabase(first, max_ring_size=2)
            with self.assertRaises(ValueError):
                lillymol_bdb.RingDatabase(first, attachment_labels=False)
            with self.assertRaises(ValueError):
                lillymol_bdb.RingDatabase([])
            with self.assertRaises(RuntimeError):
                lillymol_bdb.RingDatabase(os.path.join(tmpdir, "missing.bdb"))

    def test_import_and_construct(self):
        self.assertIsNotNone(lillymol_bdb.Selimsteg())
        precedent = lillymol_bdb.SyntheticPrecedentDatabases()
        self.assertIn("Synthetic Precedent database", repr(precedent))
        self.assertTrue(precedent.set_max_radius(2))
        precedent.close()

    def test_open_missing_database_fails(self):
        lookup = lillymol_bdb.Selimsteg()
        self.assertFalse(lookup.open_database("/no/such/berkeleydb/file.bdb"))
        lookup.close()

    def test_selimsteg_lookup(self):
        loader = _runfile_path("BerkeleyDB/iwbdb_load")
        if loader is None:
            self.skipTest("iwbdb_load not available in runfiles")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fname = os.path.join(tmpdir, "id_smiles.txt")
            dbname = os.path.join(tmpdir, "id_smiles.bdb")
            with open(input_fname, "w", encoding="utf-8") as writer:
                writer.write("mol1 CCO\n")
                writer.write("mol2 c1ccccc1\n")

            subprocess.run([loader, "-d", dbname, input_fname], check=True)

            # In a with block so the database is closed before the temporary
            # directory goes away. Leaving it open passes on a local filesystem,
            # where unlinking an open file works, but on NFS the unlink becomes a
            # silly rename to .nfsXXXX and cleanup fails with ENOTEMPTY. The bazel
            # sandbox puts /tmp on NFS here.
            with lillymol_bdb.Selimsteg() as lookup:
                self.assertTrue(lookup.open_database(dbname))
                self.assertEqual(lookup.get_smiles("mol1"), "CCO")
                self.assertIsNone(lookup.get_smiles("missing"))

                mol = lookup.get_molecule("mol2")
                self.assertIsNotNone(mol)
                self.assertEqual(mol.unique_smiles(), "c1ccccc1")
                self.assertEqual(mol.name(), "mol2")

                mols = lookup.get_molecules(["mol1", "missing", "mol2"])
                self.assertEqual(len(mols), 3)
                self.assertEqual(mols[0].unique_smiles(), "OCC")
                self.assertEqual(mols[1].natoms(), 0)
                self.assertEqual(mols[2].unique_smiles(), "c1ccccc1")

                lookup.close()

    def test_structure_database_lookup(self):
        loader = _runfile_path("Molecule_Tools_Bdb/structure_database_load")
        if loader is None:
            self.skipTest("structure_database_load not available in runfiles")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fname = os.path.join(tmpdir, "molecules.smi")
            dbstem = os.path.join(tmpdir, "structures")
            with open(input_fname, "w", encoding="utf-8") as writer:
                writer.write("CCO ethanol\n")
                writer.write("C1CCCCC1 cyclohexane\n")

            subprocess.run([loader, "-d", dbstem, input_fname], check=True)

            db = lillymol_bdb.StructureDatabase()
            self.assertTrue(db.open_read(dbstem))
            self.assertEqual(lillymol_bdb.value(lillymol_bdb.LookupParams.EXACT), 1)
            self.assertEqual(lillymol_bdb.value(lillymol_bdb.LookupParams.GRAPH), 8)

            ethanol = lillymol.MolFromSmiles("CCO query")
            self.assertEqual(db.lookup(ethanol), "ethanol")

            propanol = lillymol.MolFromSmiles("CCCO propanol")
            self.assertIsNone(db.lookup(propanol))

            strip_mask = lillymol_bdb.value(lillymol_bdb.LookupParams.STRIP)
            mixture = lillymol.MolFromSmiles("CCO.O mixture")
            self.assertEqual(db.lookup(mixture, strip_mask), "ethanol")

            graph_mask = lillymol_bdb.value(lillymol_bdb.LookupParams.GRAPH)
            self.assertEqual(db.lookup(ethanol, graph_mask), "ethanol")

            db.close()

    def test_substituent_identification_lookup(self):
        tool = _runfile_path("Molecule_Tools_Bdb/substituent_identification")
        if tool is None:
            self.skipTest("substituent_identification not available in runfiles")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fname = os.path.join(tmpdir, "build.smi")
            dbname = os.path.join(tmpdir, "substituents.bdb")
            with open(input_fname, "w", encoding="utf-8") as writer:
                writer.write("CC ethane\n")
                writer.write("CO methanol\n")
                writer.write("CN methylamine\n")

            subprocess.run(
                [tool, "-d", dbname, "-B", "-R", "2", "-w", "1", "-M", "2", input_fname],
                check=True,
            )

            with lillymol_bdb.SubstituentIdentificationLookup() as lookup:
                self.assertTrue(lookup.open_database(dbname))
                lookup.set_default_new_molecule_starting_points(1)
                lookup.set_max_substituent_size(2)

                mol = lillymol.MolFromSmiles("C methane")
                replacements = lookup.generate_replacements(mol)

                self.assertEqual(len(replacements), 3)
                self.assertEqual(
                    sorted(replacement.donor for replacement in replacements),
                    ["ethane", "methanol", "methylamine"],
                )
                for replacement in replacements:
                    self.assertEqual(replacement.name, "methane")
                    self.assertEqual(replacement.radius, 1)
                    self.assertGreaterEqual(replacement.examples, 1)
                    self.assertTrue(replacement.smiles)
                    self.assertEqual(replacement.molecule.name(), "methane")
                lookup.close()

    def test_substituent_identification_query_driven_replacement(self):
        tool = _runfile_path("Molecule_Tools_Bdb/substituent_identification")
        if tool is None:
            self.skipTest("substituent_identification not available in runfiles")

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fname = os.path.join(tmpdir, "phenyl_substituents.smi")
            dbname = os.path.join(tmpdir, "phenyl_substituents.bdb")
            with open(input_fname, "w", encoding="utf-8") as writer:
                writer.write("Clc1ccccc1 chloro\n")
                writer.write("Nc1ccccc1 nitrogen\n")
                writer.write("Oc1ccccc1 oxygen\n")
                writer.write("Sc1ccccc1 sulphur\n")

            subprocess.run(
                [tool, "-d", dbname, "-B", "-R", "4", "-Y", "dbproto", input_fname],
                check=True,
            )

            with lillymol_bdb.SubstituentIdentificationLookup() as lookup:
                self.assertTrue(lookup.open_database(dbname))
                self.assertTrue(lookup.add_query_from_smarts("[r]-!@*"))
                lookup.set_break_molecule_at_first_two_matched_atoms(True)

                phenol = lillymol.MolFromSmiles("Oc1ccccc1 phenol")
                replacements = lookup.generate_replacements(phenol)

                self.assertEqual(
                    sorted(replacement.smiles for replacement in replacements),
                    [
                        "[1C]1(=CC=CC=C1)[1Cl]",
                        "[1C]1(=CC=CC=C1)[1NH2]",
                        "[1C]1(=CC=CC=C1)[1SH]",
                    ],
                )
                self.assertEqual(
                    sorted(replacement.donor for replacement in replacements),
                    ["chloro", "nitrogen", "sulphur"],
                )
                lookup.close()

            with lillymol_bdb.SubstituentIdentificationLookup() as lookup:
                self.assertTrue(lookup.open_database(dbname))
                self.assertTrue(lookup.add_query_from_smarts("[r]-!@*"))
                lookup.set_break_molecule_at_first_two_matched_atoms(True)
                lookup.set_max_matches_per_input_molecule(2)

                phenol = lillymol.MolFromSmiles("Oc1ccccc1 phenol")
                replacements = lookup.generate_replacements(phenol)

                self.assertEqual(len(replacements), 2)
                lookup.close()

            with lillymol_bdb.SubstituentIdentificationLookup() as lookup:
                self.assertTrue(lookup.open_database(dbname))
                self.assertTrue(lookup.add_query_from_smarts("[r]-!@*"))
                lookup.set_break_molecule_at_first_two_matched_atoms(True)
                lookup.set_remove_isotopes_from_product(True)

                phenol = lillymol.MolFromSmiles("Oc1ccccc1 phenol")
                replacements = lookup.generate_replacements(phenol)

                self.assertEqual(
                    sorted(replacement.smiles for replacement in replacements),
                    [
                        "C1(=CC=CC=C1)Cl",
                        "C1(=CC=CC=C1)N",
                        "C1(=CC=CC=C1)S",
                    ],
                )
                lookup.close()

    def test_selimsteg_with_no_database_open(self):
        """Used to dereference a null Db* and take the process down."""
        lookup = lillymol_bdb.Selimsteg()
        self.assertIsNone(lookup.get_smiles("anything"))
        self.assertIsNone(lookup.get_molecule("anything"))
        self.assertEqual(lookup.get_molecules(["a", "b"])[0].natoms(), 0)
        lookup.close()          # idempotent, and safe with nothing open
        lookup.close()



if __name__ == "__main__":
    unittest.main()
