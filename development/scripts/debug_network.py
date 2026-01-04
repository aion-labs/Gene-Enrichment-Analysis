#!/usr/bin/env python3

import sys
import os
sys.path.append('code')

def test_dot_generation_direct():
    """Test DOT generation directly with mock data"""
    
    # Import the IterativeEnrichment class
    from iter_enrichment import IterativeEnrichment
    
    # Create a mock IterativeEnrichment object
    class MockIterativeEnrichment:
        def __init__(self):
            self.results = [
                {
                    "Iteration": 1,
                    "Term": "TERM1",
                    "Description": "Test term 1",
                    "Library": "Test Library",
                    "p-value": 0.01,
                    "Overlap size": "2/2",
                    "Genes": ["GENE1", "GENE2"]
                },
                {
                    "Iteration": 2,
                    "Term": "TERM2", 
                    "Description": "Test term 2",
                    "Library": "Test Library",
                    "p-value": 0.02,
                    "Overlap size": "2/2",
                    "Genes": ["GENE3", "GENE4"]
                }
            ]
        
        def to_dot(self):
            """
            Generate a valid Graphviz DOT for the iterative enrichment network,
            with sanitized, quoted IDs, semicolons, and no duplicates.
            """
            import re

            def _sanitize_id(raw: str) -> str:
                """
                Convert raw label into a valid DOT node ID: replace non-alphanumeric with underscores.
                Collapse multiple underscores and strip leading/trailing underscores.
                """
                # replace non-word characters with underscore
                s = re.sub(r"\W+", "_", raw)
                # collapse underscores
                s = re.sub(r"_+", "_", s)
                return s.strip("_")

            def _format_term_name(term_name: str) -> str:
                """
                Format term name by replacing underscores with spaces and adding colon after library shortcut.
                """
                # Import the format_term_name function from enrichment module
                from enrichment import format_term_name
                return format_term_name(term_name)

            nodes = set()
            edges = set()

            # Build nodes and edges
            for rec in self.results:
                term_label = rec.get("Term", "")
                # Format term name for display (convert underscores to spaces)
                formatted_term_label = _format_term_name(term_label)
                # sanitize and quote term ID
                raw_id = f"term_{rec['Iteration']}_{term_label}"
                term_id = _sanitize_id(raw_id)
                term_node = (
                    f'"{term_id}" '
                    f'[label="{formatted_term_label}", '
                    f'style=filled, fontcolor="white", type="term"];'
                )
                nodes.add(term_node)

                for gene in rec.get("Genes", []):
                    gene_id = _sanitize_id(f"gene_{gene}")
                    gene_node = f'"{gene_id}" [label="{gene}", type="gene"];'
                    nodes.add(gene_node)
                    # create edge with quoted IDs and semicolon
                    edge = f'"{gene_id}" -- "{term_id}";'
                    edges.add(edge)

            # Assemble DOT
            lines = []
            lines.append("graph iterative_enrichment {")
            lines.append("  graph [layout=neato];")
            lines.append("  node [shape=ellipse];")

            # add nodes
            for node in sorted(nodes):
                lines.append(f"  {node}")
            # add edges
            for edge in sorted(edges):
                lines.append(f"  {edge}")

            lines.append("}")
            return "\n".join(lines)
    
    # Test DOT generation
    try:
        mock_it = MockIterativeEnrichment()
        dot = mock_it.to_dot()
        print(f"DOT generated successfully: {len(dot)} characters")
        print("Generated DOT:")
        print(dot)
        return dot
    except Exception as e:
        print(f"Error generating DOT: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_full_network_pipeline():
    """Test the full network generation pipeline"""
    
    # Generate DOT
    dot = test_dot_generation_direct()
    if not dot:
        return False
    
    try:
        # Test merge_iterative_dot
        from code.ui.dot_utils import merge_iterative_dot
        selected_dots = {"Test Library": dot}
        merged_dot = merge_iterative_dot(selected_dots)
        print(f"Merged DOT generated: {len(merged_dot)} characters")
        
        # Test dot_to_plotly
        from code.ui.dot_utils import dot_to_plotly
        fig = dot_to_plotly(merged_dot)
        print(f"Plotly figure generated successfully: {type(fig)}")
        
        return True
    except Exception as e:
        print(f"Error in full pipeline: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Testing DOT generation...")
    test_dot_generation_direct()
    
    print("\nTesting full network pipeline...")
    test_full_network_pipeline()
