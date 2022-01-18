import os
import pytest
import sys
import json

# The scripts for this repo are in the script_dir.
# We want to import the functions we're testing from the script and add to PYTHONPATH.
# data_dir_name contains path to kraken reports needed for testing.
# We test against outputs generated by original Perl script.

this_file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir_name = os.path.join('data', 'parse-kraken-report2')
data_dir = os.path.join(this_file_dir, data_dir_name)
# path to script we are testing
script_dir = os.path.join(this_file_dir, os.pardir, os.pardir, 'bin')
sys.path.insert(1, script_dir)
sys.dont_write_bytecode = True

import parse_kraken_report2

def test_non_existing_report():
    # test non existing kraken report
    args = [0, "nonexisting.txt", "test.json", 0.5, 5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: cannot find nonexisting.txt'

def test_empty_report():
    # test empty kraken report
    args = [0, os.path.join(data_dir_name, "test_empty_report.txt"), "test.json", 0.5, 5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: %s is empty' %os.path.join(data_dir_name, "test_empty_report.txt")

def test_output_json():
    # test output file, which is not json
    args = [0, os.path.join(data_dir_name, "test_kraken_report.txt"), "test.txt", 0.5, 5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: output file test.txt must end with suffix .json'

def test_negative_num_threshold():
    # test negative num_threshold
    args = [0, os.path.join(data_dir_name, "test_kraken_report.txt"), "test.json", 0.5, -5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: -5000 is not a positive integer'

def test_negative_pct_threshold():
    # test negative pct_threshold
    pct_threshold = float(-0.5)
    args = [0, os.path.join(data_dir_name, "test_kraken_report.txt"), "test.json", pct_threshold, 5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: %f is not a positive number' %pct_threshold

def test_high_pct_threshold():
    # test >100% pct_threshold
    pct_threshold = float(105.6)
    args = [0, os.path.join(data_dir_name, "test_kraken_report.txt"), "test.json", pct_threshold, 5000]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        parse_kraken_report2.process_requirements(args)
    assert pytest_wrapped_e.type == SystemExit
    assert str(pytest_wrapped_e.value) == 'ERROR: %f is a %% and cannot be > 100' %pct_threshold

def test_read_and_parse_kraken_report():
    pct_threshold = 0.5
    num_threshold = 5000
    kraken_txt = 'test_kraken_report.txt'
    kraken_json = 'test_kraken_report.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [[142104, 10.08, 'Mycobacterium tuberculosis', '1773'], [78001, 5.53, 'Mycobacterium intracellulare', '1767']],
        [[1393167, 98.78, 'Mycobacterium', '1763']],
        [[1124355, 79.72, 'Mycobacterium tuberculosis complex', '77643'], [122223, 8.67, 'Mycobacterium avium complex (MAC)', '120793']],
        [[1404682, 99.6, 'Mycobacteriaceae', '1762']],
        2
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report2():
    pct_threshold = 0.1
    num_threshold = 1
    kraken_txt = 'test-corynebacterium.txt'
    kraken_json = 'test-corynebacterium.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [[48, 0.1, 'Human immunodeficiency virus 1', '11676']],
        [[4275, 8.9, 'Corynebacterium', '1716'], [83, 0.17, 'Mycobacterium', '1763'], [48, 0.1, 'Lentivirus', '11646']],
        [],
        [[4275, 8.9, 'Corynebacteriaceae', '1653'], [83, 0.17, 'Mycobacteriaceae', '1762'], [48, 0.1, 'Retroviridae', '11632']],
        1
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report3():
    pct_threshold = 0.5
    num_threshold = 5000
    kraken_txt = 'test-mixedmyco.txt'
    kraken_json = 'test-mixedmyco.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [[319481, 92.44, 'Mycobacterium abscessus', '36809'], [5097, 1.47, 'Mycobacterium tuberculosis', '1773']],
        [[345582, 99.99, 'Mycobacterium', '1763']],
        [[322723, 93.38, 'Mycobacterium chelonae group', '670516'], [7442, 2.15, 'Mycobacterium tuberculosis complex', '77643'], [5542, 1.6, 'Mycobacterium avium complex (MAC)', '120793']],
        [[345582, 99.99, 'Mycobacteriaceae', '1762']],
        2
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report4():
    pct_threshold = 0.5
    num_threshold = 1000
    kraken_txt = 'test-mixedmyco2.txt'
    kraken_json = 'test-mixedmyco2.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [[21302, 18.5, 'Mycobacterium chimaera', '222805'], [9316, 8.09, 'Mycobacterium tuberculosis', '1773'], [2713, 2.36, 'Mycobacterium bovis', '1765'], [1097, 0.95, 'Mycobacterium nonchromogenicum', '1782'], [1097, 0.95, 'Mycobacterium kumamotonense', '354243'], [1165, 1.01, 'Mycobacterium leprae', '1769']],
        [[114967, 99.84, 'Mycobacterium', '1763']],
        [[43595, 37.86, 'Mycobacterium avium complex (MAC)', '120793'], [13385, 11.62, 'Mycobacterium tuberculosis complex', '77643'], [2317, 2.01, 'Mycobacterium terrae complex', '1073531'], [1080, 0.94, 'Mycobacterium chelonae group', '670516']],
        [[114967, 99.84, 'Mycobacteriaceae', '1762']],
        6
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report5():
    pct_threshold = 0.5
    num_threshold = 1000
    kraken_txt = 'test-mixedmyco3.txt'
    kraken_json = 'test-mixedmyco3.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [],
        [[1456, 100.0, 'Mycobacterium', '1763']],
        [[1076, 73.9, 'Mycobacterium tuberculosis complex', '77643']],
        [[1456, 100.0, 'Mycobacteriaceae', '1762']],
        0
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report6():
    pct_threshold = 10
    num_threshold = 2000
    kraken_txt = 'test-mixedmyco4.txt'
    kraken_json = 'test-mixedmyco4.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [],
        [[7292, 80.27, 'Mycobacterium', '1763']],
        [[4514, 49.69, 'Mycobacterium tuberculosis complex', '77643']],
        [[7292, 80.27, 'Mycobacteriaceae', '1762']],
        0
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json

def test_read_and_parse_kraken_report7():
    pct_threshold = 0.5
    num_threshold = 2000
    kraken_txt = 'test-tb.txt'
    kraken_json = 'test-tb.json'

    # test reading txt kraken report
    test_kraken_report = os.path.join(data_dir_name, kraken_txt)
    got_read = parse_kraken_report2.read_kraken_report(test_kraken_report, pct_threshold, num_threshold)
    expect = (
        [[140109, 99.17, 'Mycobacterium tuberculosis', '1773']],
        [[140109, 99.17, 'Mycobacterium', '1763']],
        [[140109, 99.17, 'Mycobacterium tuberculosis complex', '77643']],
        [[140109, 99.17, 'Mycobacteriaceae', '1762']],
        1
    )
    assert got_read == expect

    # test analysing kraken report and producing json
    test_kraken_json_path = os.path.join(data_dir_name, kraken_json)
    with open(test_kraken_json_path, 'r') as f:
        test_kraken_json = json.load(f) 
    got_parse = parse_kraken_report2.parse_kraken_report(got_read[0], got_read[1], got_read[2], got_read[3], got_read[4], pct_threshold, num_threshold)
    assert got_parse == test_kraken_json