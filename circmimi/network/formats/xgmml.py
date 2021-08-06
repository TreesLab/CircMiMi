import copy
import random
from lxml import etree


class CyXGMML:
    NSMAP = {
        "dc": "http://purl.org/dc/elements/1.1/",
        "xlink": "http://www.w3.org/1999/xlink",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "cy": "http://www.cytoscape.org",
        None: "http://www.cs.rpi.edu/XGMML"
    }
    TYPE = {
        str: 'string',
        bool: 'boolean',
        list: 'list',
        int: 'integer'
    }

    def __init__(self, graph):
        self.create_root()
        self.load(graph)

    def create_root(self):
        attrs = {
            "id": str(random.randint(0, 5000)),
            "label": "circmimi_result",
            "directed": "1",
            "{http://www.cytoscape.org}documentVersion": "3.0"
        }
        self.root = etree.Element(
            'graph',
            attrib=attrs,
            nsmap=self.NSMAP
        )

    def load(self, graph):

        self.update_attrs(self.root, **graph.graph)

        # add nodes
        for node, attrs in graph.nodes(True):
            attrs = copy.deepcopy(attrs)
            id_ = attrs.pop('id', node)
            label = attrs.pop('label', node)
            data = attrs.pop('data', {})
            data_type = attrs.pop('data_type', "")
            graphics = attrs.pop('graphics', {})
            xgmml_node = self.add_node(self.root, id_, label, **attrs)

            self.add_att(xgmml_node, 'data_type', data_type)

            for k, v in data.items():
                self.add_att(xgmml_node, k, v)

            self.add_graphics(xgmml_node, **graphics)

        # add edges
        for source, target in graph.edges:
            attrs = copy.deepcopy(graph.edges[source, target])
            id_ = attrs.pop('id', "")
            label = attrs.pop('label', "{} - {}".format(source, target))
            data = attrs.pop('data', {})
            data_type = attrs.pop('data_type', "")
            graphics = attrs.pop('graphics', {})
            xgmml_edge = self.add_edge(
                self.root,
                id_,
                label,
                source,
                target,
                **attrs
            )

            self.add_att(xgmml_edge, 'data_type', data_type)
            self.add_att(xgmml_edge, 'source', source)
            self.add_att(xgmml_edge, 'target', target)

            for k, v in data.items():
                self.add_att(xgmml_edge, k, v)

            self.add_graphics(xgmml_edge, **graphics)

    def to_string(self):
        string = etree.tostring(
            self.root,
            xml_declaration=True,
            pretty_print=True,
            encoding='utf-8',
            standalone=True
        )

        return string

    @staticmethod
    def add_node(parent, id_, label, **attrs):
        node = etree.SubElement(
            parent,
            'node',
            id=id_,
            label=label,
            **attrs
        )

        return node

    @staticmethod
    def add_edge(parent, id_, label, source, target, **attrs):
        edge = etree.SubElement(
            parent,
            'edge',
            id=id_,
            label=label,
            source=source,
            target=target,
            **attrs
        )

        return edge

    @classmethod
    def add_att(cls, parent, name, value, **attrs):
        type_ = type(value)
        if type_ is int:
            value = str(value)
        elif type_ is bool:
            value = str(int(value))

        att = etree.SubElement(
            parent,
            'att',
            name=name,
            value=value,
            type=cls.TYPE[type_],
            **attrs
        )
        return att

    @staticmethod
    def add_graphics(parent, **attrs):
        graphics = etree.SubElement(parent, 'graphics', **attrs)
        return graphics

    @staticmethod
    def update_attrs(elt, **attrs):
        for k, v in attrs.items():
            elt.set(k, v)


def write_xgmml(graph, path):
    xgmml = CyXGMML(graph)

    with open(path, 'wb') as out:
        out.write(xgmml.to_string())
