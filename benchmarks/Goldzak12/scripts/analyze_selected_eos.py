"""Signed structural and representation-difference statistics."""
from collections import defaultdict
from fit_eos import aggregate


def statistics(rows, group):
    result = []
    groups = defaultdict(list)
    for row in rows:
        groups[row[group]].append(row)
    for label, values in sorted(groups.items()):
        record = {group: label, 'N': len(values)}
        for prop in ('a_A', 'B0_GPa'):
            for name, value in zip(('ME', 'MAE', 'RMSE', 'MaxAE'),
                                   aggregate([r['delta_' + prop] for r in values])):
                record[name + '_' + prop] = value
        result.append(record)
    return result
