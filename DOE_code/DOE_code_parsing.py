import ast

class CallGraphVisitor(ast.NodeVisitor):
    def __init__(self):
        self.call_graph = {}  # Stores caller-callee relationships

    def visit_FunctionDef(self, node):
        current_function = node.name
        self.call_graph[current_function] = []
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.Call):
                called_function = getattr(child.func, 'id', None) or getattr(child.func, 'attr', None)
                if called_function:
                    self.call_graph[current_function].append(called_function)
        self.generic_visit(node)

def extract_call_graph(file_path):
    with open(file_path, "r") as f:
        tree = ast.parse(f.read())
    visitor = CallGraphVisitor()
    visitor.visit(tree)
    return visitor.call_graph

# Extract the call graph from doe.py
call_graph = extract_call_graph("doe.py")
print("Call Graph:", call_graph)

from graphviz import Digraph

def visualize_call_graph(call_graph, output_file="call_graph"):
    dot = Digraph(format="png")
    for caller, callees in call_graph.items():
        dot.node(caller)
        for callee in callees:
            dot.edge(caller, callee)
    dot.render(output_file, cleanup=True)
    print(f"Call graph saved to {output_file}.png")

# Visualize the extracted call graph
visualize_call_graph(call_graph)


