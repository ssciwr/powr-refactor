import os
import subprocess
import pytest
import numpy as np

MODEL_DATA_REF = np.array(
    [
        200.48528618028934,
        170.96813898862962,
        141.27907698674588,
        113.21541313779215,
        88.20625902763122,
        67.06805146763016,
        49.991012566374565,
        36.694735175459314,
        26.638377911447826,
        19.19939157825736,
        13.787069162858446,
    ]
)


@pytest.fixture()
def run_colitest(get_chain, colitest_options):
    # need aliases and env for this to run on CI
    # run colitest
    # here we need to run with and without nonopt
    if colitest_options is None:
        colitest_options = ""
    colitest_path = "${POWR_WORK}/wrjobs/colitest1"
    colitest_command = colitest_path + colitest_options
    temp = subprocess.run(
        colitest_command,
        shell=True,
        check=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )
    print(temp.stdout)
    print(temp.stderr)
    os.system("cat ${POWR_WORK}/output/colitest1.cpr")
    yield "ran colitest"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    os.system("rm -rf ${POWR_WORK}/tmp_2day/*")
    return "Cleaned colitest tmp data"


# test the correct set-up of folder structure for jobs
def test_makechain(set_vars, get_chain):
    assert set_vars.exists()
    # now check that the scratch, output, wrdata1, wrjobs
    # dirs have been generated
    scratch_dir = set_vars / "scratch"
    output_dir = set_vars / "output"
    wrdata1_dir = set_vars / "wrdata1"
    wrjobs_dir = set_vars / "wrjobs"
    assert scratch_dir.exists()
    assert output_dir.exists()
    assert wrdata1_dir.exists()
    assert wrjobs_dir.exists()
    # now check that the correct input files have been copied into the dirs
    scratch_content = [
        "formal1",
        "modify1",
        "newdatom1",
        "newformal_cards1",
        "njn1",
        "steal1",
        "wrstart1",
        "wruniq1",
    ]
    output_content = []
    wrdata1_content = [
        "CARDS",
        "DATOM",
        "FEDAT",
        "FEDAT_FORMAL",
        "FGRID",
        "FORMAL_CARDS",
        "MODEL",
        "NEWDATOM_INPUT",
        "NEWFORMAL_CARDS_INPUT",
    ]
    wrjobs_content = [
        "coli_test",
        "formal_wrh_gen",
        "formal_wrh_xxl",
        "newformal_cards1",
        "set_repeat1",
        "tmphosts",
        "wrstart_wrh_hydro",
        "wruniq1",
        "wruniq_wrh_merged",
        "colitest1",
        "formal_wrh_hydro",
        "modify1",
        "newformal_cards_gen",
        "set_steal1",
        "wrstart1",
        "wrstart_wrh_merged",
        "wruniq_wrh_dev",
        "wruniq_wrh_vd20",
        "formal1",
        "formal_wrh_merged",
        "newdatom1",
        "njn1",
        "steal1",
        "wrstart_wrh_dev",
        "wrstart_wrh_vd20",
        "wruniq_wrh_gen",
        "wruniq_wrh_xxl",
        "formal_wrh_dev",
        "formal_wrh_vd20",
        "newdatom_gen",
        "njn_wrh_gen",
        "steal1_backup",
        "wrstart_wrh_gen",
        "wrstart_wrh_xxl",
        "wruniq_wrh_hydro",
    ]
    temp = [i.name for i in scratch_dir.iterdir()]
    assert sorted(temp) == sorted(scratch_content)
    temp = [i.name for i in output_dir.iterdir()]
    assert sorted(temp) == sorted(output_content)
    temp = [i.name for i in wrdata1_dir.iterdir()]
    assert sorted(temp) == sorted(wrdata1_content)
    temp = [i.name for i in wrjobs_dir.iterdir()]
    assert sorted(temp) == sorted(wrjobs_content)


# check that colitest run produces correct output
@pytest.mark.parametrize("colitest_options", ["", " nonopt"])
def test_colitest_run(set_vars, run_colitest):
    # check that output/colitest1.cpr is there
    # check that wrdata1/MODEL_STUDY_DONE is there
    colitest_output = set_vars / "output" / "colitest1.cpr"
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    assert colitest_output.is_file()
    assert model_output.is_file()
    # compare the output from the colitest job
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    data_np = np.fromfile(model_output, dtype=float)
    assert len(data_np) == 2246784
    # the original file length is 2233856, so if the coli job does not run
    # the test will fail and the model file length will be 2233856
    assert np.allclose(data_np[256138:256149], MODEL_DATA_REF)
