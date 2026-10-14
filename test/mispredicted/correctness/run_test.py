#!/usr/bin/env python3
"""Run isolated CLI correctness checks in the regression runner's temporary cwd."""
import os
import re
from pathlib import Path
import subprocess
import sys

exe = sys.argv[1]

def write(name, content):
    Path(name).write_text(content)

def run(args, success=True):
    result = subprocess.run([exe, *args], capture_output=True, text=True, timeout=10)
    assert result.returncode >= 0, (args, result.returncode, result.stderr)
    assert 'AddressSanitizer' not in result.stderr and 'runtime error:' not in result.stderr, result.stderr
    assert (result.returncode == 0) == success, (args, result.returncode, result.stderr)
    return result

def rows(result):
    lines = result.stdout.splitlines()
    assert lines and lines[0] == 'ID EXPT AVDIFF', result.stdout
    return {fields[0]: float(fields[2]) for fields in (line.split() for line in lines[1:])}

write('exp', 'ID measured\na 0\nb 0\nc 0\n')
write('pred1', 'ID predicted\na 2\nb 1\nc 10\n')
write('pred2', 'ID predicted\na -2\nb 1\nc 10\n')
assert rows(run(['-E', 'exp', '-d', 'pred1', 'pred2'])) == {'a': 2, 'b': 1, 'c': 10}
for count in [1, 2, 3, 4, 100]:
    for balanced in [False, True]:
        args = ['-E', 'exp', '-d', '-w', str(count)]
        if balanced:
            args += ['-b']
        result = rows(run([*args, 'pred1']))
        assert len(result) == min(count, 3), (args, result)
        assert 'c' in result, (args, result)
for predictions in ['ID predicted\na 2\nb 1\nc 10\n',
                    'ID predicted\na -2\nb -1\nc -10\n',
                    'ID predicted\na 2\nb 0\nc -10\n']:
    write('threshold_pred', predictions)
    for threshold, expected in [('0.5', {k for k, v in rows(run(['-E', 'exp', '-d', 'threshold_pred'])).items() if abs(v) >= .5}),
                                ('10', {'c'}), ('20', set())]:
        assert set(rows(run(['-E', 'exp', '-d', '-W', threshold, 'threshold_pred']))) == expected
write('one_exp', 'ID measured\na 0\n')
write('one_pred', 'ID predicted\na 2\n')
assert rows(run(['-E', 'one_exp', '-d', '-w', '1', '-b', 'one_pred'])) == {'a': -2}
write('mixed', 'ID predicted\na -10\nb -9\nc 8\n')
assert set(rows(run(['-E', 'exp', '-d', '-w', '2', '-b', 'mixed']))) == {'a', 'c'}
assert rows(run(['-E', 'exp', '-d', '-R', '^pred[12]$'])) == {'a': 2, 'b': 1, 'c': 10}
run(['-E', 'exp', '-R', '['], success=False)
run(['-E', 'exp', '-R', '^no_matching_files$'], success=False)
write('bad_pred', 'ID predicted\na 2\nb nonsense\nc 10\n')
run(['-E', 'exp', 'bad_pred'], success=False)
assert rows(run(['-E', 'exp', '-d', '-F', 'bad_pred'])) == {'a': 2, 'c': 10}
write('short_pred', 'ID predicted\na 2\nb\nc 10\n')
assert rows(run(['-E', 'exp', '-d', '-F', 'short_pred'])) == {'a': 2, 'c': 10}
write('empty_pred', 'ID predicted\n')
run(['-E', 'exp', 'empty_pred'], success=False)
write('duplicates', 'ID measured\na 0\na 100\n')
run(['-E', 'duplicates', 'one_pred'], success=False)
assert rows(run(['-E', 'duplicates', '-d', '-h', 'one_pred'])) == {'a': 2}
write('class_exp', 'ID class\na A\na B\n')
write('class_pred', 'ID class\na A\n')
assert 'a A 1 1 1' in run(['-C', '-d', '-h', '-E', 'class_exp', 'class_pred']).stdout
write('class_bad', 'ID class\na B\nb\na A\n')
assert 'a A 2 1 0.5' in run(['-C', '-d', '-h', '-F', '-E', 'class_exp', 'class_bad']).stdout
# Reject aliases before opening output; input bytes must remain intact.
original = Path('exp').read_bytes()
for output in ['exp', './exp']:
    run(['-E', 'exp', '-Y', output, 'pred1'], success=False)
    assert Path('exp').read_bytes() == original
os.link('exp', 'hardlink')
run(['-E', 'exp', '-Y', 'hardlink', 'pred1'], success=False)
assert Path('exp').read_bytes() == original
write('prediction.textproto', Path('pred1').read_text())
run(['-E', 'exp', '-T', 'prediction', 'prediction.textproto'], success=False)
assert Path('prediction.textproto').read_text() == Path('pred1').read_text()
write('files.list', 'pred1\npred2\n')
run(['-E', 'exp', '-Y', 'pred1', 'F:files.list'], success=False)
run(['-E', 'exp', '-Y', 'pred1', '-R', '^pred[12]$'], success=False)
run(['-E', 'exp', '-Y', 'same.textproto', '-T', 'same', 'pred1'], success=False)
assert not Path('same.textproto').exists()
# Reversing display must not change either filter's membership, including ties.
for selection in [[], ['-w', '2'], ['-w', '2', '-b'], ['-W', '1.5'], ['-W', '20']]:
    args = ['-E', 'exp', '-d', *selection, 'pred1', 'pred2']
    normal = run(args).stdout.splitlines()
    reverse = run(['-r', *args]).stdout.splitlines()
    assert reverse == normal[:1] + normal[:0:-1], (args, normal, reverse)
write('tie_exp', 'ID measured\nz 0\naa 0\na 0\n')
write('tie_pred', 'ID predicted\naa 1\nz 1\na 1\n')
normal = run(['-E', 'tie_exp', '-d', 'tie_pred']).stdout.splitlines()
assert [row.split()[0] for row in normal[1:]] == ['a', 'aa', 'z']
reverse = run(['-E', 'tie_exp', '-d', '-r', 'tie_pred']).stdout.splitlines()
assert reverse == normal[:1] + normal[:0:-1]
selected = run(['-E', 'tie_exp', '-d', '-w', '2', 'tie_pred']).stdout
write('tie_exp', 'ID measured\na 0\nz 0\naa 0\n')
write('tie_pred', 'ID predicted\nz 1\na 1\naa 1\n')
assert run(['-E', 'tie_exp', '-d', 'tie_pred']).stdout.splitlines() == normal
assert run(['-E', 'tie_exp', '-d', '-w', '2', 'tie_pred']).stdout == selected
# Summary and textproto ordering follow the same reversal as primary output.
for suffix, reverse_flag in [('normal', []), ('reverse', ['-r'])]:
    run(['-E', 'exp', '-d', '-Y', f'summary_{suffix}', '-T', f'proto_{suffix}',
         *reverse_flag, 'pred1', 'pred2'])
summary = Path('summary_normal').read_text().splitlines()
assert Path('summary_reverse').read_text().splitlines() == summary[:1] + summary[:0:-1]
assert Path('proto_reverse.textproto').read_text().splitlines() == Path('proto_normal.textproto').read_text().splitlines()[::-1]
write('class_tie_exp', 'ID class\nz A\na A\n')
write('class_tie_pred', 'ID class\na A\nz A\n')
normal = run(['-C', '-E', 'class_tie_exp', '-d', 'class_tie_pred']).stdout.splitlines()
assert [row.split()[0] for row in normal[1:]] == ['a', 'z']
assert run(['-C', '-E', 'class_tie_exp', '-d', '-r', 'class_tie_pred']).stdout.splitlines() == normal[:1] + normal[:0:-1]
# Opposite-sign errors must not cancel in default worst/threshold selection.
assert rows(run(['-E', 'exp', '-d', '-w', '2', 'pred1', 'pred2'])) == {'a': 2, 'c': 10}
assert rows(run(['-E', 'exp', '-d', '-W', '1.5', 'pred1', 'pred2'])) == {'a': 2, 'c': 10}
assert rows(run(['-E', 'exp', '-d', '-s', '-w', '2', 'pred1', 'pred2'])) == {'b': -1, 'c': -10}
assert rows(run(['-E', 'exp', '-d', '-s', '-W', '1.5', 'pred1', 'pred2'])) == {'c': -10}
run(['-E', 'exp', '-b', 'pred1'], success=False)
# Single-observation summaries must avoid variance(), even with assertions enabled.
for flags in [[], ['-s']]:
    run(['-E', 'one_exp', '-Y', 'single_summary', '-T', 'single_proto', *flags, 'one_pred'])
    record = Path('single_summary').read_text().splitlines()[1].split()
    assert record[0] == 'a' and record[2] == '1'
    assert float(record[-1]) == 0.0
    proto = Path('single_proto.textproto').read_text()
    assert 'npred: 1' in proto and 'std_difference: 0' in proto
# Multiple observations retain sample standard deviation of the selected metric.
write('spread1', 'ID predicted\na 2\n')
write('spread2', 'ID predicted\na -4\n')
for flags, mean, stddev in [([], 3.0, 2.0 ** .5), (['-s'], 1.0, 18.0 ** .5)]:
    run(['-E', 'one_exp', '-Y', 'spread_summary', '-T', 'spread_proto', *flags, 'spread1', 'spread2'])
    record = Path('spread_summary').read_text().splitlines()[1].split()
    assert abs(float(record[6]) - mean) < 1e-5
    assert abs(float(record[-1]) - stddev) < 1e-5
    proto = Path('spread_proto.textproto').read_text()
    value = float(re.search(r'std_difference: ([^ ]+)', proto).group(1))
    assert abs(value - stddev) < 1e-5
# Apply minimum valid-record count before selection and across all regression outputs.
write('coverage1', 'ID predicted\na 2\nb 1\nc 100\n')
write('coverage2', 'ID predicted\na 4\nb 1\n')
write('coverage3', 'ID predicted\nb 1\n')
base = ['-E', 'exp', '-d', 'coverage1', 'coverage2', 'coverage3']
assert run(['-n', '1', *base]).stdout == run(base).stdout
assert rows(run(['-n', '2', *base])) == {'a': 3, 'b': 1}
assert rows(run(['-n', '3', *base])) == {'b': 1}
assert rows(run(['-n', '4', *base])) == {}
assert rows(run(['-n', '2', '-w', '1', *base])) == {'a': 3}
assert rows(run(['-n', '2', '-W', '2', *base])) == {'a': 3}
assert set(rows(run(['-n', '2', '-w', '2', '-b', *base]))) == {'a', 'b'}
for count, eligible in [('2', {'a', 'b'}), ('4', set())]:
    run(['-n', count, '-Y', 'coverage_summary', '-T', 'coverage_proto', *base])
    records = Path('coverage_summary').read_text().splitlines()[1:]
    assert {record.split()[0] for record in records} == eligible
    assert set(re.findall(r'id: "([^"\n]+)"', Path('coverage_proto.textproto').read_text())) == eligible
write('invalid_coverage', 'ID predicted\na 2\na nonsense\nb 1\nb 1\n')
assert rows(run(['-E', 'exp', '-d', '-F', '-n', '2', 'invalid_coverage'])) == {'b': 1}
write('repeated_coverage', 'ID predicted\na 2\na 4\n')
assert rows(run(['-E', 'exp', '-d', '-n', '2', 'repeated_coverage'])) == {'a': 3}
write('coverage_class_exp', 'ID class\na A\nb B\nc C\n')
write('coverage_class1', 'ID class\na X\nb B\nc C\n')
write('coverage_class2', 'ID class\nb X\n')
class_base = ['-C', '-E', 'coverage_class_exp', '-d', 'coverage_class1', 'coverage_class2']
assert run(['-n', '1', *class_base]).stdout == run(class_base).stdout
assert run(['-n', '2', *class_base]).stdout.splitlines() == ['ID EXPT N Correct Accuracy', 'b B 2 1 0.5']
assert run(['-n', '3', *class_base]).stdout.splitlines() == ['ID EXPT N Correct Accuracy']
for invalid in ['0', '-1', 'nonsense']:
    run(['-n', invalid, *base], success=False)
print('mispredicted correctness tests passed')
