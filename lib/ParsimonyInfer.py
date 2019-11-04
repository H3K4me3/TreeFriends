import itertools as _itertools
import functools as _functools
from collections import namedtuple as _namedtuple
import numpy as _np

def _new_cost_matrix(ratio):
    levels = ["a", "c", "g", "t"]
    transversion_cost = ratio
    ans = {key:transversion_cost for key in _itertools.combinations_with_replacement(levels, 2)}
    for level in levels:
        ans[(level, level)] = 0
    transitions = [("a", "g"), ("g", "a"), ("c", "t"), ("t", "c")]
    for each in transitions:
        ans[each] = 1
    return ans

COST_MATRIX = _new_cost_matrix(1.2)

IUPAC_MAP = {
    "a": ["a"],
    "c": ["c"],
    "g": ["g"],
    "t": ["t"],
    "m": ["a", "c"],
    "r": ["a", "g"],
    "w": ["a", "t"],
    "s": ["c", "g"],
    "y": ["c", "t"],
    "k": ["g", "t"],
    "v": ["a", "c", "g"],
    "h": ["a", "c", "t"],
    "d": ["a", "g", "t"],
    "b": ["c", "g", "t"],
    "n": ["a", "c", "g", "t"], # "n" is the rare case in the database.
    "-": ["a", "c", "g", "t"],
}

EDGE_LIST = [
    ("ancestry1", "hg"        ),
    ("ancestry1", "panTro"    ),
    ("ancestry2", "ancestry1" ),
    ("ancestry2", "gorGor"    ),
    ("ancestry3", "ancestry2" ),
    ("ancestry3", "ponAbe"    ),
    ("ancestry4", "ancestry3" ),
    ("ancestry4", "rheMac"    ),
]

EdgeTuple = _namedtuple("EdgeTuple", [ "_".join(edge_a_b) for edge_a_b in EDGE_LIST])
NodeTuple = _namedtuple("NodeTuple", ['hg', 'panTro', 'gorGor', 'ponAbe', 'rheMac',
                                      'ancestry1', 'ancestry2', 'ancestry3', 'ancestry4'])

def mkNodeTuple(res):
    return NodeTuple._make((res[k] for k in NodeTuple._fields))

def specialize_ambiguous_changes(nodes):
    """
    Sample nodes argument could be
    res = {
        'chromosome': 'chr8',
        'position': 60099,
         'ref': 'A',
         'sample_number': 62784,
         'allele_number': 125568,
         'hg': 'w', 'panTro': 'a', 'gorGor': 'a', 'ponAbe': 'a', 'rheMac': 'a',
         'ancestry1': 'a', 'ancestry2': 'a', 'ancestry3': 'a', 'ancestry4': 'a'
    }
    nodes = NodeTuple._make((res[k] for k in NodeTuple._fields))
    """
    assert isinstance(nodes, NodeTuple)
    ambiguous_tree = {key: IUPAC_MAP[getattr(nodes, key)] for key in NodeTuple._fields}
    specialized_values = _itertools.product(*ambiguous_tree.values())
    specialized_keys = ambiguous_tree.keys()

    min_score = None
    maximum_parsimony_trees = None
    for values in specialized_values:
        specialized_tree = dict(zip(specialized_keys, values))
        parsimony_score = _calculate_parsimony_score(specialized_tree)
        if min_score is None or parsimony_score < min_score:
            min_score = parsimony_score
            maximum_parsimony_trees = [specialized_tree]
        elif parsimony_score == min_score:
            maximum_parsimony_trees.append(specialized_tree)
    return maximum_parsimony_trees

# Argument is a specialized tree
def _calculate_parsimony_score(tree):
    final_score = 0
    for a, b in EDGE_LIST:
        code_a, code_b = tree[a], tree[b]
        change_score = COST_MATRIX[(code_a, code_b)]
        final_score += change_score
    return final_score

# Argument is a specialized tree
def _calculate_num_changes(tree):
    s = (int(tree[a]!=tree[b]) for a, b in EDGE_LIST)
    return EdgeTuple._make(s)

def stat_edge_changes(nodes):
    trees = specialize_ambiguous_changes(nodes)
    change_tuples = [ _calculate_num_changes(tree) for tree in trees]
    return EdgeTuple._make(_np.mean(change_tuples, 0))
