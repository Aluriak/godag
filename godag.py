"""Get last gene ontology file, and save it as a predecessor>successor mapping
in json or dsv format.

usage:
    godag <output> [<idtonamefile>] [--use_id] [--json]

"""

import json
import tempfile
import itertools
import urllib.request
from collections import defaultdict

import docopt
from goatools import obo_parser


GO_OBO_URL = 'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo'


def retrieve(url:str=GO_OBO_URL) -> str:
    """Return the temporary file where retrieved file from url will be saved"""
    tmp = tempfile.NamedTemporaryFile(delete=False)
    with urllib.request.urlopen(url) as response:
        tmp.write(response.read())
    return tmp.name


def obofile_to_godag(obofile:str) -> obo_parser:
    """Get GO DAG in given obo file"""
    return obo_parser.GODag(obofile)


def godag_as_idname(obo:obo_parser) -> dict:
    """Return a mapping {go term id: go term name}"""
    idname = {}
    count_doublon = 0
    for count, (id, annot) in enumerate(obo.items(), start=1):
        if annot.id in idname:
            count_doublon += 1
        idname[annot.id] = annot.name
    print('Number of retrieved name:', len(idname),
          '(excluding #doublons: {})'.format(count_doublon))
    assert count - count_doublon == len(idname)
    return idname


def godag_as_graph(obo:obo_parser, use_id=True, *, obsoletes_key=None,
                   roots_key=None, use_special_key=False) -> dict:
    """Return a mapping {pred:succs} equivalent
    to the GO DAG found in given file.

    obsoletes_key and roots_key gives the key value used to store the obsoletes
    and roots annotation in the output dict. (default: None)
    In case of equality, obsoletes will be discarded and only roots are kept.
    If use_special_key is False, both will be discarded.

    """
    def identifier(goterm):
        return goterm.id if use_id else goterm.name
    graph = defaultdict(set)
    obsoletes, roots = set(), set()
    for count, annot in enumerate(obo.values(), start=1):
        for parent in annot.parents:
            graph[identifier(parent)].add(identifier(annot))
        if len(annot.parents) == 0:
            if annot.is_obsolete:
                obsoletes.add(identifier(annot))
            else:  # it's a root, like molecular function
                roots.add(identifier(annot))
    if use_special_key:
        graph[obsoletes_key] = frozenset(obsoletes)
        graph[roots_key] = frozenset(roots)
    print('Number of retrieved parent annotation:', len(graph))
    print('Number of obo values:', count)
    print('Total annot count:', len(set(graph.keys()) | set(itertools.chain.from_iterable(graph.values()))))
    print('Total obsoletes count:', len(obsoletes))
    print('Total roots count:', len(roots))
    if use_special_key:
        if obsoletes_key == roots_key:
            print('Only roots saved (key={})'.format(roots_key))
        else:
            print('Roots saved (key={})'.format(roots_key))
            print('Obsoletes saved (key={})'.format(Obsoletes_key))
    else:
        print('Obsoletes and roots NOT saved')
    return graph


def json_serializable_graph(graph:dict) -> dict:
    """Return a new graph, equivalent to the given one,
    where sets are replaced by lists in order to be JSON-compliant"""
    return {pred: tuple(succs) for pred, succs in graph.items()}


def graph_as_json(data:dict, filename:str, *, multiple_succs=True):
    """Save given dictionnary in given filename in JSON format

    multiple_succs -- True if a pred (key) have multiple successors

    """
    data = json_serializable_graph(data) if multiple_succs else data
    with open(filename, 'w') as fd:
        json.dump(data, fd, indent='\t', ensure_ascii=True)


def graph_as_dsv(data:dict, filename:str, *, multiple_succs=True):
    """Save given dictionnary in given filename in DSV format, as:

         pred\tsucc1\tsucc2\t…\tsuccN\n

    multiple_succs -- True if a pred (key) have multiple successors

    """
    with open(filename, 'w') as fd:
        for pred, succs in data.items():
            succs = '\t'.join(succs) if multiple_succs else succs
            fd.write('{}\t{}\n'.format(pred, succs))


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    godag = obofile_to_godag(retrieve())
    graph_to_file = graph_as_json if args['--json'] else graph_as_dsv
    # save the graph as json
    graph_to_file(godag_as_graph(godag, use_id=args['--use_id']), args['<output>'])
    print('Graph saved in ' + args['<output>'],
          ('using id' if args['--use_id'] else 'using names'))
    if args['<idtonamefile>']:
        graph_to_file(godag_as_idname(godag), args['<idtonamefile>'],
                      multiple_succs=False)
        print('Mapping id:name saved in ' + args['<idtonamefile>'])
