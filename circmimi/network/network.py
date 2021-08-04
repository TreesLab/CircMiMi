import networkx as nx
from circmimi.network.utils import read_table
from circmimi.network.formats.xgmml import write_xgmml


class CyNetwork:
    def __init__(self):
        self.graph = nx.DiGraph()

    def load_data(self, file_, k1, k2, k3):
        for data in read_table(file_):
            circRNA = data[k1 - 1]
            mediator = data[k2 - 1]
            target_gene = data[k3 - 1]

            self.graph.add_node(
                circRNA,
                data_type="circRNA"
            )

            self.graph.add_node(
                mediator,
                data_type="mediator"
            )

            self.graph.add_node(
                target_gene,
                data_type='target_gene'
            )

            self.graph.add_edge(
                circRNA,
                mediator,
                data_type="circRNA-mediator"
            )

            self.graph.add_edge(
                mediator,
                target_gene,
                data_type="mediator-target_gene"
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
