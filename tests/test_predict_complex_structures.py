from __future__ import annotations

import copy
import importlib.util
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SKILL_ROOT = ROOT / "skills" / "predict-complex-structures"
SCRIPTS_ROOT = SKILL_ROOT / "scripts"
PREPARE_PATH = SCRIPTS_ROOT / "prepare_input.py"
RUNNER_PATH = SCRIPTS_ROOT / "run_protenix.py"
VALIDATOR_PATH = SCRIPTS_ROOT / "validate_output.py"
ASSET_BUILDER_PATH = SCRIPTS_ROOT / "build_runtime_assets.py"
FIXED_TIME = "2026-08-12T00:00:00+00:00"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


PREPARE = load_module("predict_complex_prepare", PREPARE_PATH)
VALIDATOR = load_module("predict_complex_validator", VALIDATOR_PATH)
RUNNER = load_module("predict_complex_runner", RUNNER_PATH)
ASSET_BUILDER = load_module("predict_complex_asset_builder", ASSET_BUILDER_PATH)


def add_fingerprint(document, key):
    document[key] = PREPARE.document_fingerprint(document)
    return document


def structure_inchikey(smiles):
    toolkit = PREPARE.load_rdkit()
    molecule = toolkit["Chem"].MolFromSmiles(smiles)
    assert molecule is not None
    return toolkit["Chem"].MolToInchiKey(molecule)


def identity_artifact(
    *,
    resolution_id="ligand-1",
    disposition="ready_for_standardization",
    handoff_status="ready",
    structure="CCO",
):
    inchikey = structure_inchikey(structure)
    standardizer_path = (
        "/" + "Users" + "/example/workspace/standardize_structures.py"
    )
    document = {
        "schema_version": "1.0.0",
        "workflow": "chemical-identity-resolution",
        "generated_at_utc": FIXED_TIME,
        "options": {"standardizer_script": standardizer_path},
        "resolutions": [
            {
                "request": {"id": resolution_id, "query": "ethanol"},
                "record_alignment_status": (
                    "exact" if handoff_status == "ready" else "ambiguous"
                ),
                "sample_identity_status": "not_assessed",
                "disposition": disposition,
                "candidates": [{"candidate_id": "candidate-1"}],
                "standardization_handoff": {
                    "status": handoff_status,
                    "target_skill": "standardize-chemical-structures",
                    "records": (
                        [
                            {
                                "id": resolution_id,
                                "structure": structure,
                                "source_candidate_id": "candidate-1",
                                "source_inchikey": inchikey,
                            }
                        ]
                        if handoff_status == "ready"
                        else []
                    ),
                },
            }
        ],
    }
    return add_fingerprint(document, "result_fingerprint")


def standardization_artifact(
    *,
    record_id="ligand-1",
    disposition="ready_for_downstream",
    standardized_structure="CCO",
    original_structure=None,
    parent_structure=None,
    review_reasons=None,
    parse_status="success",
    standardization_status="completed",
):
    original = original_structure or standardized_structure or "CCO"
    parent = parent_structure or standardized_structure
    inchikey = structure_inchikey(standardized_structure or original)
    parent_inchikey = structure_inchikey(parent) if parent else None
    document = {
        "schema_version": "1.0.0",
        "workflow": "chemical-structure-standardization-qc",
        "generated_at_utc": FIXED_TIME,
        "tool_versions": {"rdkit": "2025.9.2"},
        "options": {"profile": "chembl-pipeline"},
        "records": [
            {
                "id": record_id,
                "original_structure": original,
                "standardized_structure": standardized_structure,
                "parent_structure": parent,
                "inchikey": inchikey,
                "parent_inchikey": parent_inchikey,
                "parse_status": parse_status,
                "standardization_status": standardization_status,
                "disposition": disposition,
                "human_review_required": list(review_reasons or []),
            }
        ],
        "duplicate_groups": [],
    }
    return add_fingerprint(document, "result_fingerprint")


def base_request(*, proteins=None, ligand=None, evidence_artifact=None):
    value = {
        "schema_version": "1.0.0",
        "job_id": "complex-001",
        "research_question": "生成可复核的结构假设",
        "proteins": proteins
        or [
            {
                "chain_id": "A",
                "sequence": "ACDEFGHIK",
                "organism": "Homo sapiens",
                "source": {
                    "database": "UniProt",
                    "identifier": "P12345",
                    "url": "https://www.uniprot.org/uniprotkb/P12345",
                },
                "modifications": [],
            }
        ],
        "model": {
            "name": "protenix-v2",
            "version": "2.0.0",
            "weights_id": "protenix-v2",
            "seeds": [101],
            "msa_strategy": "none",
            "template_strategy": "none",
            "rna_msa_strategy": "none",
            "dtype": "bf16",
            "samples": 1,
            "steps": 20,
            "cycles": 4,
        },
        "covalent_bonds": [],
        "constraint": None,
    }
    if ligand is not None:
        value["ligand"] = ligand
    if evidence_artifact is not None:
        value["evidence_artifact"] = evidence_artifact
    return value


def ligand_request(identity_path, standardization_path):
    return {
        "chain_id": "L",
        "reported_name": "ethanol",
        "source": {"database": "PubChem", "identifier": "702"},
        "identity_artifact": str(identity_path),
        "identity_resolution_id": "ligand-1",
        "standardization_artifact": str(standardization_path),
        "standardization_record_id": "ligand-1",
        "input_mode": "standardized_smiles",
    }


def prepare_request(root, request):
    result_dir = root / "result"
    result, exit_code = PREPARE.prepare_job(
        request,
        root,
        result_dir,
        generated_at_utc=FIXED_TIME,
    )
    return result_dir, result, exit_code


def write_upstream_artifacts(root, identity=None, standardization=None):
    identity_path = root / "identity.json"
    standardization_path = root / "standardization.json"
    identity_path.write_text(
        json.dumps(identity or identity_artifact()),
        encoding="utf-8",
    )
    standardization_path.write_text(
        json.dumps(standardization or standardization_artifact()),
        encoding="utf-8",
    )
    return identity_path, standardization_path


def write_mmcif(path, chains):
    rows = []
    atom_id = 1
    entity_id = 1
    x = 0.0
    for chain_id, residues in chains.items():
        seq_id = 1
        for residue_name, atoms, group in residues:
            for atom_name, element in atoms:
                rows.append(
                    " ".join(
                        [
                            group,
                            str(atom_id),
                            element,
                            atom_name,
                            ".",
                            residue_name,
                            chain_id,
                            str(entity_id),
                            str(seq_id),
                            "?",
                            f"{x:.3f}",
                            "0.000",
                            "0.000",
                            "1.00",
                            "20.00",
                            "?",
                            str(seq_id),
                            residue_name,
                            chain_id,
                            atom_name,
                            "1",
                        ]
                    )
                )
                atom_id += 1
                x += 1.0
            seq_id += 1
        entity_id += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
        + "\n".join(rows)
        + "\n#\n",
        encoding="utf-8",
    )


def valid_confidence(chain_count=2):
    pair_iptm = [
        [0.0 if row == column else 0.72 for column in range(chain_count)]
        for row in range(chain_count)
    ]
    return {
        "ranking_score": 0.8,
        "plddt": 88.0,
        "ptm": 0.75,
        "iptm": 0.72,
        "gpde": 2.5,
        "chain_ptm": [0.75] * chain_count,
        "chain_iptm": [0.72] * chain_count,
        "chain_plddt": [88.0] * chain_count,
        "chain_pair_iptm": pair_iptm,
        "chain_pair_plddt": [[88.0] * chain_count for _ in range(chain_count)],
        "has_clash": False,
        "disorder": 0.0,
        "num_recycles": 4,
    }


def expected_protein_residues(count=9):
    names = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "GLY",
        "HIS",
        "ILE",
        "LYS",
    ]
    return [
        (name, [("N", "N"), ("CA", "C")], "ATOM")
        for name in names[:count]
    ]


def write_runtime_assets(root, *, model_name="protenix-v2", version="2.0.0"):
    runtime_root = root / "protenix-runtime"
    checkpoint = runtime_root / "checkpoint" / f"{model_name}.pt"
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_bytes(b"test-checkpoint")
    cache_paths = {
        "ccd_components_file": "common/components.cif",
        "ccd_components_rdkit_mol_file": "common/components.cif.rdkit_mol.pkl",
        "obsolete_release_data_csv": "common/obsolete_release_date.csv",
        "pdb_cluster_file": "common/clusters-by-entity-40.txt",
    }
    caches = []
    for name, relative in cache_paths.items():
        path = runtime_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"test cache: {name}\n", encoding="utf-8")
        caches.append(
            {
                "name": name,
                "relative_path": relative,
                "sha256": RUNNER.sha256_file(path),
            }
        )
    manifest = {
        "schema_version": "1.0.0",
        "model_name": model_name,
        "protenix_version": version,
        "protenix_root_dir": str(runtime_root),
        "checkpoint": {
            "name": "model_checkpoint",
            "relative_path": checkpoint.relative_to(runtime_root).as_posix(),
            "sha256": RUNNER.sha256_file(checkpoint),
        },
        "required_caches": caches,
    }
    manifest_path = root / "runtime-assets.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return manifest_path


def write_fake_protenix(root):
    executable = root / "protenix"
    executable.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        "  echo 'protenix, version 2.0.0'\n"
        "  exit 0\n"
        "fi\n"
        "exit 0\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def add_completed_outputs(result_dir, *, chains=None, confidence=None):
    prediction_dir = (
        result_dir / "protenix_output" / "input" / "complex-001" / "seed_101" / "predictions"
    )
    structure_path = prediction_dir / "complex-001_sample_0.cif"
    write_mmcif(
        structure_path,
        chains
        or {
            "A": expected_protein_residues(),
            "L": [
                (
                    "ETO",
                    [("C1", "C"), ("C2", "C"), ("O1", "O")],
                    "HETATM",
                )
            ],
        },
    )
    confidence_path = prediction_dir / "complex-001_summary_confidence_sample_0.json"
    confidence_path.write_text(
        json.dumps(confidence or valid_confidence()),
        encoding="utf-8",
    )
    manifest_path = result_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "status": "partial",
            "exit_code": 0,
            "failure": None,
            "outputs": {
                "structure_files": [
                    structure_path.relative_to(result_dir).as_posix()
                ],
                "auxiliary_structure_files": [],
                "summary_confidence_files": [
                    confidence_path.relative_to(result_dir).as_posix()
                ],
                "full_confidence_files": [],
            },
        }
    )
    manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return structure_path, confidence_path


class PrepareInputTests(unittest.TestCase):
    def test_skill_metadata_and_required_files(self):
        self.assertTrue((SKILL_ROOT / "SKILL.md").is_file())
        self.assertTrue((SKILL_ROOT / "agents" / "openai.yaml").is_file())
        for name in ("prepare_input.py", "run_protenix.py", "validate_output.py"):
            self.assertTrue((SCRIPTS_ROOT / name).is_file())
        text = (SKILL_ROOT / "SKILL.md").read_text(encoding="utf-8")
        self.assertIn("name: \"predict-complex-structures\"", text)
        self.assertIn("dataPro_search", text)
        self.assertNotIn("X-Agent-Plan-Key:", text)

    def test_result_directory_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            actual = root / "actual-result"
            actual.mkdir()
            linked = root / "linked-result"
            linked.symlink_to(actual, target_is_directory=True)
            with self.assertRaisesRegex(PREPARE.InputFailure, "不允许 symlink"):
                PREPARE.prepare_job(
                    base_request(),
                    root,
                    linked,
                    generated_at_utc=FIXED_TIME,
                )

    def test_input_artifact_paths_cannot_escape_request_directory(self):
        for path_factory in (
            lambda root, outside: str(outside / "evidence.json"),
            lambda root, outside: "../outside/evidence.json",
        ):
            with tempfile.TemporaryDirectory() as directory:
                workspace = Path(directory)
                root = workspace / "request"
                root.mkdir()
                outside = workspace / "outside"
                outside.mkdir()
                (outside / "evidence.json").write_text(
                    json.dumps(
                        {
                            "code": 0,
                            "dataset_type": "academic_search",
                            "items": [],
                        }
                    ),
                    encoding="utf-8",
                )
                request = base_request(
                    evidence_artifact=path_factory(root, outside)
                )
                with self.assertRaisesRegex(
                    PREPARE.InputFailure,
                    "必须位于 request 文件目录内",
                ):
                    prepare_request(root, request)

    def test_precomputed_files_require_hash_and_controlled_extension(self):
        for suffix, include_hash, expected_error in (
            (".a3m", False, "sha256 必须提供"),
            (".txt", True, "文件扩展名必须属于"),
        ):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                msa_path = root / f"protein{suffix}"
                msa_path.write_text(
                    ">query\nACDEFGHIK\n",
                    encoding="utf-8",
                )
                descriptor = {"path": str(msa_path)}
                if include_hash:
                    descriptor["sha256"] = PREPARE.sha256_file(msa_path)
                request = base_request()
                request["proteins"][0]["unpaired_msa"] = descriptor
                request["model"]["msa_strategy"] = "precomputed"
                with self.assertRaisesRegex(
                    PREPARE.InputFailure,
                    expected_error,
                ):
                    prepare_request(root, request)

    def test_prepare_protein_only_official_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, result, exit_code = prepare_request(root, base_request())
            self.assertEqual(exit_code, 0)
            task = result["protenix_input"][0]
            self.assertEqual(task["name"], "complex-001")
            self.assertEqual(
                task["sequences"][0],
                {
                    "proteinChain": {
                        "sequence": "ACDEFGHIK",
                        "count": 1,
                        "id": ["A"],
                        "modifications": [],
                    }
                },
            )
            self.assertEqual(
                result["ligand_handoff"]["ligand_status"],
                "not_applicable",
            )
            self.assertTrue((result_dir / "input_manifest.json").is_file())

    def test_prepare_protein_ligand_uses_upstream_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            identity_path, standardization_path = write_upstream_artifacts(root)
            request = base_request(
                ligand=ligand_request(identity_path, standardization_path)
            )
            _, result, exit_code = prepare_request(root, request)
            self.assertEqual(exit_code, 0)
            ligand = result["protenix_input"][0]["sequences"][1]["ligand"]
            self.assertEqual(ligand["ligand"], "CCO")
            handoff = result["ligand_handoff"]
            self.assertEqual(handoff["selected_form"]["selection"], "standardized_form")
            self.assertFalse(
                handoff["selected_form"]["parent_used_for_prediction"]
            )
            self.assertEqual(
                handoff["expected_output_profile"]["heavy_atom_count"],
                3,
            )
            self.assertEqual(
                handoff["expected_output_profile"]["formal_charge"],
                0,
            )
            self.assertIn(
                "upstream_review_reasons",
                handoff["selected_form"],
            )
            self.assertEqual(
                handoff["upstream"]["identity"]["handoff_status"],
                "ready",
            )
            self.assertEqual(handoff["artifact_binding"]["status"], "matched")
            serialized = json.dumps(result, ensure_ascii=False)
            self.assertNotIn("/" + "Users" + "/example", serialized)

    def test_prepare_multi_chain_input(self):
        proteins = [
            {
                "chain_id": "A",
                "sequence": "ACDE",
                "organism": "human",
                "source": {"database": "UniProt", "identifier": "P1"},
            },
            {
                "chain_id": "B",
                "sequence": "FGHIK",
                "organism": "human",
                "source": {"database": "UniProt", "identifier": "P2"},
            },
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, result, _ = prepare_request(
                root,
                base_request(proteins=proteins),
            )
            ids = [
                item["proteinChain"]["id"][0]
                for item in result["protenix_input"][0]["sequences"]
            ]
            self.assertEqual(ids, ["A", "B"])

    def test_prepare_dna_and_rna_require_provenance(self):
        request = base_request()
        request["dna"] = [
            {
                "chain_id": "D",
                "sequence": "ATGC",
                "organism": "synthetic construct",
                "source": {"database": "user-provided", "identifier": "dna-1"},
            }
        ]
        request["rna"] = [
            {
                "chain_id": "R",
                "sequence": "AUGC",
                "organism": "synthetic construct",
                "source": {"database": "user-provided", "identifier": "rna-1"},
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            _, result, exit_code = prepare_request(Path(directory), request)
            self.assertEqual(exit_code, 0)
            self.assertEqual(
                [item["chain_id"] for item in result["input_manifest"]["polymers"]],
                ["A", "D", "R"],
            )
            self.assertEqual(
                result["input_manifest"]["polymers"][1]["entity_type"],
                "dnaSequence",
            )

    def test_dna_without_source_fails_closed(self):
        request = base_request()
        request["dna"] = [
            {
                "chain_id": "D",
                "sequence": "ATGC",
                "organism": "synthetic construct",
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, "source.database"):
                prepare_request(Path(directory), request)

    def test_invalid_protein_sequence_fails_closed(self):
        request = base_request()
        request["proteins"][0]["sequence"] = "ACDZ"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(PREPARE.InputFailure, "非法字符"):
                prepare_request(root, request)

    def test_uncontrolled_model_with_extra_weights_is_rejected(self):
        request = base_request()
        request["model"]["name"] = "protenix_mini_esm_v0.5.0"
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                PREPARE.InputFailure,
                "额外权重模型尚未支持",
            ):
                prepare_request(Path(directory), request)

    def test_template_strategy_rejects_incompatible_model(self):
        request = base_request()
        request["model"]["name"] = "protenix_base_default_v0.5.0"
        request["model"]["template_strategy"] = "precomputed"
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                PREPARE.InputFailure,
                "不支持模板推理",
            ):
                prepare_request(Path(directory), request)

    def test_invalid_modification_position_fails_closed(self):
        request = base_request()
        request["proteins"][0]["modifications"] = [
            {"ptmType": "CCD_MSE", "ptmPosition": 99}
        ]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, "超出实体范围"):
                prepare_request(Path(directory), request)

    def test_unbounded_model_workload_fails_closed(self):
        for key, value in (
            ("samples", 101),
            ("steps", 10_001),
            ("cycles", 101),
        ):
            request = base_request()
            request["model"][key] = value
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(PREPARE.InputFailure, "必须是"):
                    prepare_request(Path(directory), request)

        request = base_request()
        request["model"]["seeds"] = [2**32]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, r"2\^32"):
                prepare_request(Path(directory), request)

    def test_invalid_covalent_bond_entity_fails_closed(self):
        request = base_request()
        request["covalent_bonds"] = [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": "1",
                "atom1": "CA",
                "entity2": 2,
                "copy2": 1,
                "position2": "1",
                "atom2": "C1",
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, "超出实体范围"):
                prepare_request(Path(directory), request)

    def test_weights_id_must_match_actual_model_selection(self):
        request = base_request()
        request["model"]["weights_id"] = "protenix_base_default_v1.0.0"
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, "weights_id"):
                prepare_request(Path(directory), request)

    def test_covalent_bond_forces_scientific_review(self):
        request = base_request()
        request["covalent_bonds"] = [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": "1",
                "atom1": "CA",
                "entity2": 1,
                "copy2": 1,
                "position2": "2",
                "atom2": "CA",
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            _, result, exit_code = prepare_request(Path(directory), request)
            self.assertEqual(exit_code, 0)
            self.assertEqual(
                result["run_manifest"]["review_disposition"],
                "review_required",
            )
            self.assertIn(
                "COVALENT_BOND_NOT_INDEPENDENTLY_VERIFIED",
                result["input_manifest"]["scientific_review_reasons"],
            )

    def test_valid_constraint_is_normalized_for_constraint_model(self):
        request = base_request()
        request["model"]["name"] = "protenix_base_constraint_v0.5.0"
        request["model"]["weights_id"] = "protenix_base_constraint_v0.5.0"
        request["constraint"] = {
            "contact": [
                {
                    "entity1": 1,
                    "copy1": 1,
                    "position1": 2,
                    "atom1": "CA",
                    "entity2": 1,
                    "copy2": 1,
                    "position2": 8,
                    "max_distance": 8,
                    "min_distance": 3,
                }
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            _, result, exit_code = prepare_request(Path(directory), request)
            self.assertEqual(exit_code, 0)
            contact = result["protenix_input"][0]["constraint"]["contact"][0]
            self.assertEqual(contact["max_distance"], 8.0)
            self.assertEqual(contact["min_distance"], 3.0)
            self.assertEqual(
                result["run_manifest"]["review_disposition"],
                "review_required",
            )
            self.assertIn(
                "SOFT_CONSTRAINT_NOT_INDEPENDENTLY_VERIFIED",
                result["input_manifest"]["scientific_review_reasons"],
            )

    def test_invalid_constraint_distance_fails_closed(self):
        request = base_request()
        request["model"]["name"] = "protenix_base_constraint_v0.5.0"
        request["model"]["weights_id"] = "protenix_base_constraint_v0.5.0"
        request["constraint"] = {
            "contact": [
                {
                    "entity1": 1,
                    "copy1": 1,
                    "position1": 2,
                    "entity2": 1,
                    "copy2": 1,
                    "position2": 8,
                    "max_distance": 3,
                    "min_distance": 8,
                }
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(PREPARE.InputFailure, "min_distance"):
                prepare_request(Path(directory), request)

    def test_unconfirmed_identity_blocks_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            identity_path, standardization_path = write_upstream_artifacts(
                root,
                identity=identity_artifact(
                    disposition="review_required",
                    handoff_status="blocked_pending_resolution",
                ),
            )
            result_dir, result, exit_code = prepare_request(
                root,
                base_request(
                    ligand=ligand_request(identity_path, standardization_path)
                ),
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(result["protenix_input"], [])
            self.assertFalse(result["run_manifest"]["inference_allowed"])
            self.assertEqual(
                result["ligand_handoff"]["review_disposition"],
                "review_required",
            )
            self.assertTrue((result_dir / "ligand_handoff.json").is_file())

    def test_reported_ligand_without_identity_artifact_blocks_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, standardization_path = write_upstream_artifacts(root)
            ligand = ligand_request(root / "missing.json", standardization_path)
            ligand.pop("identity_artifact")
            _, result, exit_code = prepare_request(
                root,
                base_request(ligand=ligand),
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(result["protenix_input"], [])
            self.assertIn(
                "R-LIGAND-IDENTITY-ARTIFACT-MISSING",
                result["ligand_handoff"]["review_reasons"],
            )

    def test_precomputed_msa_requires_a_file_for_each_protein(self):
        request = base_request()
        request["model"]["msa_strategy"] = "precomputed"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(
                PREPARE.InputFailure,
                "每条蛋白链提供",
            ):
                prepare_request(root, request)

    def test_secret_in_precomputed_msa_is_rejected_before_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            msa_path = root / "protein.a3m"
            msa_path.write_text(
                ">query\nACDEFGHIK\napi_key=secret-review-value\n",
                encoding="utf-8",
            )
            request = base_request()
            request["proteins"][0]["unpaired_msa"] = {
                "path": str(msa_path),
                "sha256": PREPARE.sha256_file(msa_path),
            }
            request["model"]["msa_strategy"] = "precomputed"
            with self.assertRaisesRegex(
                PREPARE.InputFailure,
                "possible credential",
            ):
                prepare_request(root, request)

    def test_identity_and_standardization_mismatch_blocks_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            aspirin = "CC(=O)Oc1ccccc1C(=O)O"
            identity_path, standardization_path = write_upstream_artifacts(
                root,
                identity=identity_artifact(structure="CCO"),
                standardization=standardization_artifact(
                    standardized_structure=aspirin,
                    original_structure=aspirin,
                ),
            )
            _, result, exit_code = prepare_request(
                root,
                base_request(
                    ligand=ligand_request(identity_path, standardization_path)
                ),
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(result["protenix_input"], [])
            self.assertEqual(
                result["ligand_handoff"]["artifact_binding"]["status"],
                "mismatch",
            )
            self.assertIn(
                "R-UPSTREAM-ARTIFACT-MISMATCH",
                result["ligand_handoff"]["review_reasons"],
            )

    def test_stereochemistry_review_propagates_without_using_parent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lactic_acid = "C[C@H](O)C(=O)O"
            identity_path, standardization_path = write_upstream_artifacts(
                root,
                identity=identity_artifact(structure=lactic_acid),
                standardization=standardization_artifact(
                    disposition="review_required",
                    standardized_structure=lactic_acid,
                    original_structure=lactic_acid,
                    review_reasons=["R-STEREO-CONFLICT"],
                ),
            )
            _, result, exit_code = prepare_request(
                root,
                base_request(
                    ligand=ligand_request(identity_path, standardization_path)
                ),
            )
            self.assertEqual(exit_code, 0)
            handoff = result["ligand_handoff"]
            self.assertEqual(handoff["review_disposition"], "review_required")
            self.assertIn("R-STEREO-CONFLICT", handoff["review_reasons"])
            self.assertIn(
                "@",
                result["protenix_input"][0]["sequences"][1]["ligand"]["ligand"],
            )
            self.assertFalse(
                handoff["selected_form"]["parent_used_for_prediction"]
            )

    def test_rejected_standardization_blocks_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            identity_path, standardization_path = write_upstream_artifacts(
                root,
                standardization=standardization_artifact(
                    disposition="rejected",
                    standardized_structure=None,
                    parse_status="error",
                    standardization_status="not_run",
                ),
            )
            _, result, exit_code = prepare_request(
                root,
                base_request(
                    ligand=ligand_request(identity_path, standardization_path)
                ),
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(
                result["ligand_handoff"]["review_disposition"],
                "rejected",
            )
            self.assertEqual(result["protenix_input"], [])

    def test_datapro_not_configured_is_not_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, result, _ = prepare_request(root, base_request())
            self.assertEqual(result["evidence"]["evidence_status"], "not_run")
            self.assertEqual(result["evidence"]["items"], [])

    def test_datapro_missing_doi_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence_path = root / "evidence.json"
            evidence = {
                "code": 0,
                "msg": "success",
                "trace_id": "trace-public-test",
                "dataset_type": "academic_search",
                "query": "protein ligand",
                "items": [
                    {
                        "name": "Paper without DOI",
                        "url": "https://example.org/paper",
                        "abstract": "Candidate evidence only.",
                        "authors": "A. Author",
                        "date_published": "2026-01-01",
                        "doi": "",
                    }
                ],
            }
            evidence_path.write_text(json.dumps(evidence), encoding="utf-8")
            _, result, _ = prepare_request(
                root,
                base_request(evidence_artifact=str(evidence_path)),
            )
            self.assertEqual(result["evidence"]["evidence_status"], "completed")
            self.assertFalse(result["evidence"]["item_review"][0]["doi_present"])
            self.assertEqual(result["evidence"]["response"]["items"], evidence["items"])

    def test_datapro_nested_doi_is_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence_path = root / "evidence.json"
            evidence = {
                "code": 0,
                "msg": "success",
                "trace_id": "trace-public-shape",
                "dataset_type": "academic_search",
                "items": [
                    {
                        "name": "Protenix technical report",
                        "url": "https://doi.org/10.1101/2025.01.08.631967",
                        "extra_data": {
                            "doi": "https://doi.org/10.1101/2025.01.08.631967"
                        },
                    }
                ],
            }
            evidence_path.write_text(json.dumps(evidence), encoding="utf-8")
            result_dir, result, exit_code = prepare_request(
                root,
                base_request(evidence_artifact=str(evidence_path)),
            )
            self.assertEqual(exit_code, 0)
            self.assertEqual(result["evidence"]["evidence_status"], "completed")
            self.assertTrue(result["evidence"]["item_review"][0]["doi_present"])
            self.assertEqual(
                result["evidence"]["response"]["items"], evidence["items"]
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(
                any("缺少 DOI" in item for item in report["warnings"])
            )

    def test_repeated_preparation_has_stable_fingerprints(self):
        first_request = base_request()
        second_request = copy.deepcopy(first_request)
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            _, first, _ = prepare_request(Path(first_dir), first_request)
            _, second, _ = prepare_request(Path(second_dir), second_request)
            self.assertEqual(
                first["input_manifest"]["input_fingerprint"],
                second["input_manifest"]["input_fingerprint"],
            )
            self.assertEqual(
                first["run_manifest"]["run_fingerprint"],
                second["run_manifest"]["run_fingerprint"],
            )


class RunnerTests(unittest.TestCase):
    def setUp(self):
        RUNNER._CANCEL_REQUESTED = False

    def test_official_command_uses_manifest_parameters(self):
        manifest = {
            "model": {
                "name": "protenix-v2",
                "seeds": [101, 102],
                "msa_strategy": "none",
                "template_strategy": "precomputed",
                "rna_msa_strategy": "none",
                "dtype": "bf16",
                "samples": 2,
                "steps": 20,
                "cycles": 4,
            }
        }
        command = RUNNER.build_command("/runtime/protenix", manifest)
        self.assertEqual(command[:6], ["/runtime/protenix", "pred", "-i", "protenix_input.json", "-o", "protenix_output"])
        self.assertIn("101,102", command)
        self.assertEqual(command[command.index("--use_template") + 1], "true")
        self.assertEqual(command[command.index("--use_msa") + 1], "false")

    def test_official_wounresol_structure_is_classified_as_auxiliary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main = (
                root
                / "protenix_output"
                / "input"
                / "complex-001"
                / "seed_101"
                / "predictions"
                / "complex-001_sample_0.cif"
            )
            write_mmcif(main, {"A": expected_protein_residues()})
            auxiliary = main.with_name("complex-001_sample_0_wounresol.cif")
            auxiliary.write_text(main.read_text(encoding="utf-8"), encoding="utf-8")
            outputs = RUNNER.discover_outputs(root)
            self.assertEqual(len(outputs["structure_files"]), 1)
            self.assertEqual(len(outputs["auxiliary_structure_files"]), 1)

    def test_protenix_not_installed_is_not_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            report, exit_code = RUNNER.run_job(
                result_dir,
                protenix_command="definitely-not-installed-protenix",
                generated_at_utc=FIXED_TIME,
            )
            manifest = json.loads(
                (result_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(manifest["status"], "not_run")
            self.assertEqual(
                manifest["failure"]["code"],
                "PROTENIX_NOT_INSTALLED",
            )
            self.assertTrue(report["valid"], report)

    def test_gpu_unavailable_is_not_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            fake = root / "protenix"
            fake.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
            with patch.object(RUNNER, "detect_gpu", return_value=None):
                report, exit_code = RUNNER.run_job(
                    result_dir,
                    protenix_command=str(fake),
                    generated_at_utc=FIXED_TIME,
                )
            manifest = json.loads(
                (result_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(manifest["failure"]["code"], "GPU_UNAVAILABLE")
            self.assertTrue(report["valid"], report)

    def test_runtime_assets_manifest_is_required_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            fake = write_fake_protenix(root)
            with patch.object(
                RUNNER,
                "detect_gpu",
                return_value={"devices": [{"name": "test-gpu"}]},
            ):
                report, exit_code = RUNNER.run_job(
                    result_dir,
                    protenix_command=str(fake),
                    generated_at_utc=FIXED_TIME,
                )
            manifest = json.loads(
                (result_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(exit_code, 2)
            self.assertEqual(
                manifest["failure"]["code"],
                "ASSETS_MANIFEST_REQUIRED",
            )
            self.assertTrue(report["valid"], report)

    def test_runner_rejects_tampered_protenix_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            input_path = result_dir / "protenix_input.json"
            document = json.loads(input_path.read_text(encoding="utf-8"))
            document[0]["sequences"][0]["proteinChain"]["sequence"] = "AAAAAAAAA"
            input_path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(
                RUNNER.RunFailure,
                "protenix_input.json 与 input_manifest 未绑定",
            ):
                RUNNER.run_job(
                    result_dir,
                    protenix_command="definitely-not-installed-protenix",
                    generated_at_utc=FIXED_TIME,
                )

    def test_runner_rejects_tampered_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            evidence_path = result_dir / "evidence.json"
            document = json.loads(evidence_path.read_text(encoding="utf-8"))
            document["notice"] = "tampered"
            document["evidence_fingerprint"] = PREPARE.document_fingerprint(
                document
            )
            evidence_path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(
                RUNNER.RunFailure,
                "evidence.json 与 input_manifest 未绑定",
            ):
                RUNNER.run_job(
                    result_dir,
                    protenix_command="definitely-not-installed-protenix",
                    generated_at_utc=FIXED_TIME,
                )

    def test_runner_rejects_model_drift_from_prepared_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            manifest_path = result_dir / "run_manifest.json"
            document = json.loads(manifest_path.read_text(encoding="utf-8"))
            document["model"]["steps"] = 200
            document["run_fingerprint"] = PREPARE.document_fingerprint(document)
            manifest_path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(
                RUNNER.RunFailure,
                "run_manifest.model 与 input_manifest 未绑定",
            ):
                RUNNER.run_job(
                    result_dir,
                    protenix_command="definitely-not-installed-protenix",
                    generated_at_utc=FIXED_TIME,
                )

    def test_runtime_assets_manifest_validates_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = write_runtime_assets(root)
            summary, runtime_root = RUNNER.validate_runtime_assets(
                manifest_path,
                base_request()["model"],
            )
            self.assertEqual(summary["model_name"], "protenix-v2")
            self.assertFalse(summary["downloads_allowed"])
            self.assertEqual(
                summary["checkpoint"]["sha256"],
                RUNNER.sha256_file(
                    runtime_root / summary["checkpoint"]["relative_path"]
                ),
            )
            self.assertNotIn(str(runtime_root), json.dumps(summary))

    def test_runtime_assets_manifest_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = write_runtime_assets(root)
            linked = root / "runtime-assets-link.json"
            linked.symlink_to(manifest_path)
            with self.assertRaisesRegex(
                RUNNER.RunFailure,
                "缺失或不是普通文件",
            ):
                RUNNER.validate_runtime_assets(
                    linked,
                    base_request()["model"],
                )

    def test_runner_log_symlink_is_rejected_without_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "outside.txt"
            target.write_text("keep", encoding="utf-8")
            (root / "protenix_stdout.txt").symlink_to(target)
            with self.assertRaisesRegex(RUNNER.RunFailure, "不允许 symlink"):
                RUNNER.write_logs(root, "new", "")
            self.assertEqual(target.read_text(encoding="utf-8"), "keep")

    def test_runtime_assets_builder_matches_runner_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_manifest = write_runtime_assets(root)
            source = json.loads(source_manifest.read_text(encoding="utf-8"))
            runtime_root = Path(source["protenix_root_dir"])
            built = ASSET_BUILDER.build_manifest(
                runtime_root,
                "protenix-v2",
                "2.0.0",
            )
            output = root / "built-runtime-assets.json"
            ASSET_BUILDER.write_manifest(output, built)
            summary, validated_root = RUNNER.validate_runtime_assets(
                output,
                base_request()["model"],
            )
            self.assertEqual(validated_root, runtime_root.resolve())
            self.assertEqual(summary["protenix_version"], "2.0.0")

    def test_runtime_assets_builder_fails_on_missing_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_manifest = write_runtime_assets(root)
            source = json.loads(source_manifest.read_text(encoding="utf-8"))
            runtime_root = Path(source["protenix_root_dir"])
            (runtime_root / "common" / "components.cif").unlink()
            with self.assertRaisesRegex(ASSET_BUILDER.AssetFailure, "文件缺失"):
                ASSET_BUILDER.build_manifest(
                    runtime_root,
                    "protenix-v2",
                    "2.0.0",
                )

    def test_runtime_assets_builder_rejects_symlinked_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_manifest = write_runtime_assets(root)
            source = json.loads(source_manifest.read_text(encoding="utf-8"))
            runtime_root = Path(source["protenix_root_dir"])
            checkpoint = runtime_root / "checkpoint" / "protenix-v2.pt"
            target = root / "outside-checkpoint.pt"
            target.write_bytes(checkpoint.read_bytes())
            checkpoint.unlink()
            checkpoint.symlink_to(target)
            with self.assertRaisesRegex(ASSET_BUILDER.AssetFailure, "symlink"):
                ASSET_BUILDER.build_manifest(
                    runtime_root,
                    "protenix-v2",
                    "2.0.0",
                )

    def test_runner_success_requires_validated_assets_and_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            fake = write_fake_protenix(root)
            assets = write_runtime_assets(root)

            def fake_execute(
                _command,
                cwd,
                _timeout,
                _cancel_file,
                environment,
            ):
                self.assertEqual(
                    Path(environment["PROTENIX_ROOT_DIR"]).resolve(),
                    (root / "protenix-runtime").resolve(),
                )
                self.assertNotIn("AGENT_PLAN_API_KEY", environment)
                prediction_dir = (
                    cwd
                    / "protenix_output"
                    / "input"
                    / "complex-001"
                    / "seed_101"
                    / "predictions"
                )
                structure_path = prediction_dir / "complex-001_sample_0.cif"
                write_mmcif(
                    structure_path,
                    {
                        "A": [
                            *expected_protein_residues()
                        ]
                    },
                )
                confidence_path = (
                    prediction_dir
                    / "complex-001_summary_confidence_sample_0.json"
                )
                confidence_path.write_text(
                    json.dumps(valid_confidence(chain_count=1)),
                    encoding="utf-8",
                )
                return {
                    "exit_code": 0,
                    "stdout": "completed",
                    "stderr": "",
                    "timed_out": False,
                    "canceled": False,
                    "elapsed_seconds": 1.0,
                }

            with (
                patch.object(
                    RUNNER,
                    "detect_gpu",
                    return_value={"devices": [{"name": "test-gpu"}]},
                ),
                patch.object(RUNNER, "execute_process", side_effect=fake_execute),
                patch.dict(
                    os.environ,
                    {"AGENT_PLAN_API_KEY": "secret-not-forwarded"},
                ),
            ):
                report, exit_code = RUNNER.run_job(
                    result_dir,
                    protenix_command=str(fake),
                    runtime_assets_manifest=assets,
                    generated_at_utc=FIXED_TIME,
                )
            self.assertEqual(exit_code, 0)
            self.assertTrue(report["completion_ready"], report)
            manifest = json.loads(
                (result_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "completed")
            self.assertFalse(
                manifest["runtime"]["assets"]["downloads_allowed"]
            )

    def test_runner_rejects_reused_directory_with_stale_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result_dir, _, _ = prepare_request(root, base_request())
            stale = result_dir / "protenix_output" / "stale.cif"
            stale.parent.mkdir(parents=True)
            stale.write_text("stale", encoding="utf-8")
            fake = write_fake_protenix(root)
            with self.assertRaisesRegex(RUNNER.RunFailure, "陈旧运行产物"):
                RUNNER.run_job(
                    result_dir,
                    protenix_command=str(fake),
                    generated_at_utc=FIXED_TIME,
                )

    def test_process_timeout_is_reported(self):
        result = RUNNER.execute_process(
            [sys.executable, "-c", "import time; time.sleep(2)"],
            Path.cwd(),
            timeout_seconds=0.1,
        )
        self.assertTrue(result["timed_out"])
        self.assertFalse(result["canceled"])

    def test_cancel_file_stops_process(self):
        with tempfile.TemporaryDirectory() as directory:
            cancel_file = Path(directory) / "cancel"
            cancel_file.touch()
            result = RUNNER.execute_process(
                [sys.executable, "-c", "import time; time.sleep(2)"],
                Path.cwd(),
                timeout_seconds=5,
                cancel_file=cancel_file,
            )
            self.assertTrue(result["canceled"])
            self.assertFalse(result["timed_out"])

    def test_nonzero_process_exit_is_preserved(self):
        result = RUNNER.execute_process(
            [sys.executable, "-c", "raise SystemExit(7)"],
            Path.cwd(),
            timeout_seconds=5,
        )
        self.assertEqual(result["exit_code"], 7)
        self.assertFalse(result["timed_out"])


class ValidatorTests(unittest.TestCase):
    def prepared_complex(self, root):
        identity_path, standardization_path = write_upstream_artifacts(root)
        return prepare_request(
            root,
            base_request(
                ligand=ligand_request(identity_path, standardization_path)
            ),
        )[0]

    def test_valid_completed_fixture_is_ready_for_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(result_dir)
            report = VALIDATOR.validate_result_dir(
                result_dir,
                generated_at_utc=FIXED_TIME,
            )
            self.assertTrue(report["valid"], report)
            self.assertTrue(report["completion_ready"])
            self.assertEqual(
                report["structure_checks"][0]["ligand_check"]["status"],
                "matched",
            )

    def test_official_auxiliary_wounresol_file_does_not_break_pairing(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            structure_path, _ = add_completed_outputs(result_dir)
            auxiliary = structure_path.with_name(
                "complex-001_sample_0_wounresol.cif"
            )
            auxiliary.write_text(
                structure_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["auxiliary_structure_files"] = [
                auxiliary.relative_to(result_dir).as_posix()
            ]
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(report["completion_ready"], report)

    def test_confirmed_multicomponent_ccd_requires_exact_residue_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            request = base_request(
                ligand={
                    "chain_id": "L",
                    "reported_name": "test glycan",
                    "source": {"database": "PDB CCD", "identifier": "NAG_BMA"},
                    "input_mode": "ccd",
                    "ccd_id": "CCD_NAG_BMA",
                    "confirmation_status": "expert_confirmed",
                }
            )
            result_dir, _, _ = prepare_request(root, request)
            add_completed_outputs(
                result_dir,
                chains={
                    "A": expected_protein_residues(),
                    "L": [
                        ("NAG", [("C1", "C")], "HETATM"),
                        ("BMA", [("C1", "C")], "HETATM"),
                    ],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(report["completion_ready"], report)

    def test_missing_mmcif_blocks_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            _, confidence_path = add_completed_outputs(result_dir)
            for path in result_dir.rglob("*.cif"):
                path.unlink()
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["structure_files"] = []
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            self.assertTrue(confidence_path.is_file())
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("缺少 CIF/mmCIF" in item for item in report["errors"])
            )

    def test_missing_confidence_blocks_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(result_dir)
            for path in result_dir.rglob("*summary_confidence*.json"):
                path.unlink()
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["summary_confidence_files"] = []
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("summary confidence" in item for item in report["errors"])
            )

    def test_invalid_mmcif_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            structure_path, _ = add_completed_outputs(result_dir)
            structure_path.write_text("not mmcif", encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["valid"])
            self.assertTrue(
                any("mmCIF 无法解析" in item for item in report["errors"])
            )

    def test_chain_mapping_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(
                result_dir,
                chains={
                    "B": [("ALA", [("N", "N"), ("CA", "C")], "ATOM")],
                    "L": [
                        (
                            "ETO",
                            [("C1", "C"), ("C2", "C"), ("O1", "O")],
                            "HETATM",
                        )
                    ],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["valid"])
            self.assertTrue(any("缺少期望链" in item for item in report["errors"]))

    def test_unexpected_chain_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(
                result_dir,
                chains={
                    "A": [
                        *expected_protein_residues()
                    ],
                    "L": [
                        (
                            "ETO",
                            [("C1", "C"), ("C2", "C"), ("O1", "O")],
                            "HETATM",
                        )
                    ],
                    "X": [("ALA", [("N", "N"), ("CA", "C")], "ATOM")],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("包含未请求链" in item for item in report["errors"])
            )

    def test_truncated_polymer_chain_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(
                result_dir,
                chains={
                    "A": [
                        *expected_protein_residues(8)
                    ],
                    "L": [
                        (
                            "ETO",
                            [("C1", "C"), ("C2", "C"), ("O1", "O")],
                            "HETATM",
                        )
                    ],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("链 A 残基数不一致" in item for item in report["errors"])
            )

    def test_same_length_wrong_polymer_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            wrong_sequence = expected_protein_residues()
            wrong_sequence[3] = (
                "ALA",
                [("N", "N"), ("CA", "C")],
                "ATOM",
            )
            add_completed_outputs(
                result_dir,
                chains={
                    "A": wrong_sequence,
                    "L": [
                        (
                            "ETO",
                            [("C1", "C"), ("C2", "C"), ("O1", "O")],
                            "HETATM",
                        )
                    ],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any(
                    "链 A 残基序列与输入不一致" in item
                    for item in report["errors"]
                )
            )

    def test_declared_protein_modification_is_validated_as_modified_residue(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            request = base_request()
            request["proteins"][0]["modifications"] = [
                {"ptmType": "CCD_MSE", "ptmPosition": 2}
            ]
            result_dir, _, _ = prepare_request(root, request)
            modified = expected_protein_residues()
            modified[1] = (
                "MET",
                [("N", "N"), ("CA", "C")],
                "ATOM",
            )
            add_completed_outputs(
                result_dir,
                chains={"A": modified},
                confidence={
                    "ranking_score": 0.8,
                    "plddt": 88.0,
                    "ptm": 0.75,
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(report["completion_ready"], report)

    def test_unlisted_prediction_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            structure_path, _ = add_completed_outputs(result_dir)
            hidden = structure_path.with_name("complex-001_sample_1.cif")
            hidden.write_text(
                structure_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("与磁盘实际输出不一致" in item for item in report["errors"])
            )

    def test_stale_generated_report_does_not_poison_revalidation(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(result_dir)
            (result_dir / "validation_report.json").write_text(
                json.dumps({"api_key": "stale-secret-value"}),
                encoding="utf-8",
            )
            (result_dir / "report.html").write_text(
                "/" + "Users" + "/stale/private/path",
                encoding="utf-8",
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(report["completion_ready"], report)

    def test_invalid_full_confidence_is_rejected_when_present(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            structure_path, _ = add_completed_outputs(result_dir)
            full_path = structure_path.with_name(
                "complex-001_full_data_sample_0.json"
            )
            full_path.write_text('{"atom_plddt": NaN}', encoding="utf-8")
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["full_confidence_files"] = [
                full_path.relative_to(result_dir).as_posix()
            ]
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("含非有限数值" in item for item in report["errors"])
            )

    def test_out_of_range_and_wrong_shape_confidence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            confidence = valid_confidence()
            confidence["plddt"] = 101.0
            confidence["ptm"] = "high"
            confidence["chain_pair_iptm"] = [[0.0]]
            add_completed_outputs(result_dir, confidence=confidence)
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("plddt 高于 100.0" in item for item in report["errors"])
            )
            self.assertTrue(
                any("ptm 必须是有限数值" in item for item in report["errors"])
            )
            self.assertTrue(
                any("chain_pair_iptm 必须是 2x2 矩阵" in item for item in report["errors"])
            )

    def test_optional_official_confidence_fields_are_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            confidence = valid_confidence()
            confidence["has_clash"] = "false"
            confidence["num_recycles"] = -1
            confidence["chain_ptm"] = [0.7]
            confidence["chain_pair_plddt"][0][1] = 101.0
            add_completed_outputs(result_dir, confidence=confidence)
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("has_clash 必须是 boolean" in item for item in report["errors"])
            )
            self.assertTrue(
                any("num_recycles 必须是非负整数" in item for item in report["errors"])
            )
            self.assertTrue(
                any("chain_ptm 必须是长度 2" in item for item in report["errors"])
            )
            self.assertTrue(
                any("chain_pair_plddt" in item for item in report["errors"])
            )

    def test_scientific_review_reason_cannot_be_removed_from_risky_input(self):
        request = base_request()
        request["covalent_bonds"] = [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": "1",
                "atom1": "CA",
                "entity2": 1,
                "copy2": 1,
                "position2": "2",
                "atom2": "CA",
            }
        ]
        with tempfile.TemporaryDirectory() as directory:
            result_dir, _, _ = prepare_request(Path(directory), request)
            for name, fingerprint_key in (
                ("input_manifest.json", "input_fingerprint"),
                ("run_manifest.json", "run_fingerprint"),
            ):
                path = result_dir / name
                document = json.loads(path.read_text(encoding="utf-8"))
                document["scientific_review_reasons"] = []
                document[fingerprint_key] = PREPARE.document_fingerprint(document)
                path.write_text(json.dumps(document), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["valid"])
            self.assertTrue(
                any("与共价键/软约束输入不一致" in item for item in report["errors"])
            )

    def test_ligand_atom_identity_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            add_completed_outputs(
                result_dir,
                chains={
                    "A": [
                        *expected_protein_residues()
                    ],
                    "L": [("ETO", [("C1", "C"), ("O1", "O")], "HETATM")],
                },
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["valid"])
            self.assertTrue(
                any("配体重原子数不一致" in item for item in report["errors"])
            )

    def test_secret_and_sensitive_path_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            sensitive_path = "/" + "Users" + "/example/private/file"
            (result_dir / "protenix_stdout.txt").write_text(
                "Bearer " + "A" * 32 + "\n" + sensitive_path,
                encoding="utf-8",
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["valid"])
            self.assertTrue(
                any("credential" in item for item in report["errors"])
            )
            self.assertTrue(
                any("host-sensitive path" in item for item in report["errors"])
            )
            (result_dir / "protenix_stdout.txt").write_text(
                '{"api_key":"test-secret-value"}',
                encoding="utf-8",
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(
                any("credential" in item for item in report["errors"])
            )
            (result_dir / "protenix_stdout.txt").write_text(
                "ark-" + "a" * 48 + "\n"
                "/mlx_devbox/users/example/project/result.json",
                encoding="utf-8",
            )
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertTrue(
                any("credential" in item for item in report["errors"])
            )
            self.assertTrue(
                any("host-sensitive path" in item for item in report["errors"])
            )

    def test_structure_and_confidence_counts_must_match(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            structure_path, _ = add_completed_outputs(result_dir)
            extra_structure = structure_path.with_name("complex-001_sample_1.cif")
            extra_structure.write_text(
                structure_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["structure_files"].append(
                extra_structure.relative_to(result_dir).as_posix()
            )
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertTrue(
                any("文件数量不一致" in item for item in report["errors"])
            )

    def test_structure_and_confidence_sample_keys_must_match(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            _, confidence_path = add_completed_outputs(result_dir)
            renamed = confidence_path.with_name(
                "complex-001_summary_confidence_sample_1.json"
            )
            confidence_path.rename(renamed)
            manifest_path = result_dir / "run_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["outputs"]["summary_confidence_files"] = [
                renamed.relative_to(result_dir).as_posix()
            ]
            manifest["run_fingerprint"] = PREPARE.document_fingerprint(manifest)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            report = VALIDATOR.validate_result_dir(result_dir)
            self.assertFalse(report["completion_ready"])
            self.assertEqual(report["output_pairing"]["status"], "error")
            self.assertTrue(
                any("job/seed/sample 不一致" in item for item in report["errors"])
            )

    def test_upstream_identity_and_standardization_fingerprints_are_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            identity_path, standardization_path = write_upstream_artifacts(root)
            result_dir, result, _ = prepare_request(
                root,
                base_request(
                    ligand=ligand_request(identity_path, standardization_path)
                ),
            )
            handoff = result["ligand_handoff"]
            self.assertEqual(
                handoff["upstream"]["identity"]["result_fingerprint"],
                identity_artifact()["result_fingerprint"],
            )
            self.assertEqual(
                handoff["upstream"]["standardization"]["result_fingerprint"],
                standardization_artifact()["result_fingerprint"],
            )
            report = VALIDATOR.validate_result_dir(
                result_dir,
                generated_at_utc=FIXED_TIME,
            )
            self.assertTrue(report["valid"], report)

    def test_validator_cli_runs_independently(self):
        with tempfile.TemporaryDirectory() as directory:
            result_dir = self.prepared_complex(Path(directory))
            completed = subprocess.run(
                [
                    sys.executable,
                    str(VALIDATOR_PATH),
                    str(result_dir),
                    "--generated-at",
                    FIXED_TIME,
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((result_dir / "validation_report.json").is_file())
            self.assertTrue((result_dir / "report.html").is_file())


if __name__ == "__main__":
    unittest.main()
