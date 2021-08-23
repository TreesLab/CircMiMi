import networkx as nx
from circmimi.network.utils import read_table
from circmimi.network.formats.xgmml import write_xgmml


class CyNetwork:
    def __init__(self):
        self.graph = nx.DiGraph()

    def load_data(self, file_, k1, k2, k3, header=True):
        for data in read_table(file_, header=header):
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

    def update_node_attr(self, node, key, val):
        self.graph.nodes[node][key] = val

    def update_edge_attr(self, edge, key, val):
        self.graph.edges[edge][key] = val

    def apply_layout(self, layout):
        layout.apply(self.graph)

    def apply_style(self, style):
        style.apply(self.graph)

    def to_xgmml(self, file_):
        write_xgmml(self.graph, file_)


class Layout:
    def __init__(self):
        pass

    def apply(self, graph):
        pass


class Style:
    def __init__(self):
        self.circRNAs_style = {'type': 'circle', 'fill': '#F292AF'}
        self.mediator_style = {'type': 'ROUNDED_RECTANGLE', 'fill': '#8BE26E'}
        self.target_gene_style = {'type': 'diamond', 'fill': '#6D9AF4'}

        self.circRNAs_mediator_style = {'fill': '#FFC1C1'}
        self.mediator_taget_gene_style = {'fill': '#9BDDA2'}

    def apply(self, graph):
        for node in graph.nodes:
            node_type = graph.nodes[node]['data_type']
            if node_type == 'circRNA':
                graph.nodes[node]['graphics'] = self.circRNAs_style
            elif node_type == 'mediator':
                graph.nodes[node]['graphics'] = self.mediator_style
            elif node_type == 'target_gene':
                graph.nodes[node]['graphics'] = self.target_gene_style

        for edge in graph.edges:
            edge_type = graph.edges[edge]['data_type']
            if edge_type == 'circRNA-mediator':
                graph.edges[edge]['graphics'] = self.circRNAs_mediator_style
            elif edge_type == 'mediator-target_gene':
                graph.edges[edge]['graphics'] = self.mediator_taget_gene_style
