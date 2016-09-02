"""Get last gene ontology file, and save it as a predecessor>successor mapping
in json format

usage:
    godag <output> [<idtonamefile>] [--use_id]

"""

import json
import tempfile
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
    for id, annot in obo.items():
        idname[annot.id] = annot.name
    return idname


def godag_as_graph(obo:obo_parser, use_id=True) -> dict:
    """Return a mapping {pred:succs} equivalent
    to the GO DAG found in given file"""
    def identifier(goterm):
        return goterm.id if use_id else goterm.name
    graph = defaultdict(set)
    for annot in obo.values():
        for parent in annot.parents:
            graph[identifier(parent)].add(identifier(annot))
    return graph


def serializable_graph(graph:dict) -> dict:
    """Return a new graph, equivalent to the given one,
    where sets are replaced by lists"""
    return {pred: tuple(succs) for pred, succs in graph.items()}


def graph_as_json(data:dict, filename:str):
    """Save given dictionnary in given filename in JSON format"""
    with open(filename, 'w') as fd:
        json.dump(data, fd)


if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    godag = obofile_to_godag(retrieve())
    # save the graph as json
    graph_as_json(serializable_graph(godag_as_graph(godag, use_id=args['--use_id'])), args['<output>'])
    print('Graph saved in ' + args['<output>'],
          ('using id' if args['--use_id'] else 'using names'))
    if args['<idtonamefile>']:
        graph_as_json(godag_as_idname(godag), args['<idtonamefile>'])
        print('Mapping id:name saved in ' + args['<idtonamefile>'])
