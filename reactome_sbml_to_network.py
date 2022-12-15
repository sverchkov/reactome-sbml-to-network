"""A python tool for extracting network structure from the SBML files provided by Reactome."""

from pathlib import Path
from xml.etree import ElementTree
import logging

import pandas as pd

LOG = logging.getLogger()

# Assuming that namespaces used in the SBML files stay consistent
NS = {
    'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
    'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'bqbiol': 'http://biomodels.net/biology-qualifiers/'
}

ID_PREFIX = 'https://identifiers.org/'


def parse_collection(dir_path, nodes_csv = None, edges_csv = None):
    """Parses a collection (a directory of SBML files) into nodes and edges dataframes and optionally saves them to CSVs"""

    # Parse individual files
    dfs = [parse_sbml(fpath) for fpath in Path(dir_path).iterdir()]

    # Unzip results
    node_dfs, edge_dfs = zip(*dfs)

    # Concatenate and remove duplicates
    nodes_df = pd.concat(node_dfs, ignore_index=True).drop_duplicates(ignore_index=True)
    edges_df = pd.concat(edge_dfs, ignore_index=True).drop_duplicates(ignore_index=True)

    # Optionally save
    if nodes_csv: nodes_df.to_csv(nodes_csv, index=None)
    if edges_csv: edges_df.to_csv(edges_csv, index=None)

    return nodes_df, edges_df


def parse_sbml(file_path):
    """Parse a single SBML file into a pair of dataframes (nodes and edges)"""

    root = ElementTree.parse(file_path).getroot()

    # Parse reactions
    reactions = root.findall('./sbml:model/sbml:listOfReactions/sbml:reaction', NS)

    reactions_parsed = {
        reaction.attrib['id']: {
            'reactants': get_species_ids(reaction.find('./sbml:listOfReactants', NS)),
            'products': get_species_ids(reaction.find('./sbml:listOfProducts', NS)),
            'modifiers': get_species_ids(reaction.find('./sbml:listOfModifiers', NS))
        }
        for reaction in reactions
    }

    # Characterize species
    species_mapping = {}
    complexes = {}

    species_set = {x for reaction in reactions_parsed.values() for group in reaction.values() for x in group}

    for species_id in species_set:
        species_element = root.find(f"./sbml:model/sbml:listOfSpecies/sbml:species[@id='{species_id}']", NS)
        if not species_element:
            LOG.warning(f'Could not find {species_id} in species list')
        else:
            the_is_element = species_element.find('./sbml:annotation/rdf:RDF/rdf:Description/bqbiol:is', NS)
            reactome_id = get_rdf_bag_resource_ids(the_is_element, 'reactome')[0]
            uniprot_ids = get_rdf_bag_resource_ids(the_is_element, 'uniprot')
            chebi_ids = get_rdf_bag_resource_ids(the_is_element, 'chebi')
            has_part_element = species_element.find('./sbml:annotation/rdf:RDF/rdf:Description/bqbiol:hasPart', NS)
            if has_part_element: # This is a complex
                species_mapping[species_id] = {'id': reactome_id, 'type': 'complex'}
                complexes[reactome_id] = get_rdf_bag_resource_ids(has_part_element, 'uniprot')
            elif len(uniprot_ids) == 1: # This is a protein
                species_mapping[species_id] = {'id': uniprot_ids[0], 'type': 'protein'}
            elif len(chebi_ids) == 1: # This is a small molecule
                species_mapping[species_id] = {'id': chebi_ids[0], 'type': 'molecule'}
            else: # Unknown ID
                species_mapping[species_id] = {'id': reactome_id, 'type': 'unknown'}
    
    
    # Table for physical entities
    entities_df = pd.DataFrame.from_dict(species_mapping, orient='index')

    # Table for reactions
    reactions_df = pd.DataFrame({'id': reactions_parsed.keys()})
    reactions_df['pathway'] = get_rdf_bag_resource_ids(
        root.find('./sbml:model/sbml:annotation/rdf:RDF/rdf:Description/bqbiol:is', NS),
        'reactome'
    )[0]
    reactions_df['type'] = 'reaction'

    # Table for reaction edges
    reaction_edges = pd.DataFrame.from_records(
        [
            {'entity': entity, 'relation': relation, 'target': reaction_id}
            for reaction_id, reaction in reactions_parsed.items()
            for relation, entities in reaction.items()
            for entity in entities
        ]
    )

    # Remap entity IDs
    try:
        reaction_edges['source'] = reaction_edges.entity.map(entities_df.id)
    except:
        print(entities_df)
        print(entities_df.index.duplicated())
        print(reaction_edges)
        print(reaction_edges.index.duplicated())
        raise
    
    # Rename relations
    reaction_edges['relation'] = reaction_edges.relation.map({
        'reactants': 'reactant_of',
        'products': 'product_of',
        'modifiers': 'modifier_of'
    })

    # Table for complex edges
    complex_edges = pd.DataFrame.from_records([
        {'source': protein, 'relation': 'part_of', 'target': complex_id}
        for complex_id, complex_members in complexes.items() for protein in complex_members
    ])

    # Compose node table
    nodes_df = pd.concat([entities_df, reactions_df], ignore_index=True)

    # Compose edge table
    edges_df = pd.concat([reaction_edges[['source', 'relation', 'target']], complex_edges], ignore_index=True)

    return nodes_df, edges_df


def get_species_ids(species_ref_element_container):
    """Gets the species ids from the reaction members lists."""

    if species_ref_element_container:
        return [
            element.attrib['species']
            for element in species_ref_element_container.findall('./sbml:speciesReference', NS)
        ]
    else:
        return []


def strip_prefix(in_str, prefix):
    """Strip a prefix from a string only if it is an exact match."""

    if in_str[0:len(prefix)] == prefix:
        return in_str[len(prefix):]
    else:
        LOG.warning(f'Could not strip prefix {prefix} from string {in_str}')
        return in_str


def get_rdf_bag_resource_ids(member_element, resource_name):
    """Get the list of resource IDs from the RDF bag in the indicated element"""

    rsc = f"{{{NS['rdf']}}}resource"
    return [
        strip_prefix(entry.attrib[rsc], ID_PREFIX)
        for entry in member_element.findall('./rdf:Bag/rdf:li', NS)
        if entry.attrib[rsc].startswith(ID_PREFIX+resource_name+':')
    ]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='A python tool for extracting network structure from the SBML files provided by Reactome.'
    )

    parser.add_argument('sbml_dir')
    parser.add_argument('nodes_csv')
    parser.add_argument('edges_csv')

    args = parser.parse_args()

    parse_collection(args.sbml_dir, args.nodes_csv, args.edges_csv)