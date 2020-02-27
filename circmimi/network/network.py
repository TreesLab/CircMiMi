import networkx as nx
from circmimi.network.utils import read_table
from circmimi.network.formats.xgmml import write_xgmml


class CyNetwork:
    def __init__(self):
        self.graph = nx.DiGraph()

    def load_data(self, file_):
        for data in read_table(file_):
            circRNA_id = ':'.join(data[:5])

            self.graph.add_node(
                circRNA_id,
                data_type="circRNA",
                data=dict(list(data._asdict().items())[:5])
            )

            self.graph.add_node(
                data.mirna,
                data_type="miRNA",
                data={'name': data.mirna}
            )

            self.graph.add_node(
                data.target_gene,
                data_type='target_gene',
                data={'name': data.target_gene}
            )

            self.graph.add_edge(
                circRNA_id,
                data.mirna,
                data_type="circ-mir",
                data=dict(list(data._asdict().items())[6:9])
            )

            self.graph.add_edge(
                data.mirna,
                data.target_gene,
                data_type="mir-target",
                data=dict(list(data._asdict().items())[10:])
            )

    def update_node(self, node, key_val):
        self.graph.nodes[node].update(key_val)

    def update_edge(self, edge, key_val):
        self.graph.edges[edge].update(key_val)

    def apply_layout(self, layout):
        pass

    def apply_style(self, style):
        pass

    def to_xgmml(self, file_):
        write_xgmml(self.graph, file_)


class Layout:
    def __init__(self):
        pass

    def apply(self, graph):
        return graph


class Style:
    def __init__(self):
        pass

    def apply(self, graph):
        return graph
