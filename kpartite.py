# %%
import os
import re
import warnings
import pandas as pd
import numpy as np
import rpy2.robjects as robjects

def find_data(dir, apply_threshold=None, default_threshold=0.001, debug=False):
    """ Find all the csv files in the data/ directory. If applying a threshold, we search for any files with that substring
    If files contain multiple substrings from the apply_threshold dictionary, the last one will be used.

    Args:
        dir (string): The directory to search
        apply_threshold (dictionary, optional): An optional mapping if certain files should have different thresholds. Defaults to None.

    Returns:
        list: A list of dicts containing dataframes and their metadata

        Note: Had to change from just a list of dataframes because
        R doesn't understand arbitrary user-defined dataframe attributes
    """

    dataframe_datalist = []
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            if debug: print(f"Loaded {file}")
            data = {}
            path = os.path.join(dir, file)
            df = pd.read_csv(path, index_col=0)
            data['name'] = file
            data['anti'] = True if 'anti' in file.lower() else False
            data['threshold'] = default_threshold
            data['dataframe'] = df
            if apply_threshold is not None:
                for key, value in apply_threshold.items():
                    if key in file:
                        data['threshold'] = value
            if debug: print(f"Threshold for {file} is {data['threshold']}\n")
            dataframe_datalist.append(data)
    return dataframe_datalist

def convert_mapped(mapped_results, apply_threshold=None, default_threshold=0.001, debug=False):

    # Create an empty list to store dictionaries
    dataframe_datalist = []

    # Iterate over each element in the original list
    for file, double in mapped_results.items():
        print('yay1\n\n\n')

        # Convert numpy array to pandas DataFrame
        df = pd.DataFrame(double)

        # Access R row and column names using reticulate
        r_matrix = robjects.r['matrix'](double)
        df.columns = robjects.r.colnames(r_matrix)
        df.index = robjects.r.rownames(r_matrix)

        # Append the DataFrame to the list
        dataframe_datalist.append(df)

        # Create a dictionary with 'name' and 'dataframe' keys
        data = {'name': file,
                    'dataframe': df,
                    'anti': True if 'anti' in file.lower() else False,
                    'threshold': default_threshold
                    }
        if apply_threshold is not None:
            for key, value in apply_threshold.items():
                    if key in file:
                        data['threshold'] = value
        if debug: print(f"Threshold for {file} is {data['threshold']}\n")
        if debug: print(f'Columns are {data["dataframe"].columns}')
        # Append the dictionary to the list_of_dictionaries
        dataframe_datalist.append(data)

    return dataframe_datalist


# %%
def filter_df(data, selected=None, value=None, anti=None, drop='higher', keep='lower', debug=False):
    """"Filter a dataframe by a particular column and thresholding value; return an edgelist

    Args:
        df (Pandas dataframe): The dataframe of (anti/)associations
        selected (string, optional): The particular column of interest; leave blank or set to None if not filtering
        value (float, optional): The value to threshold by
        anti (bool, optional): Whether it is an anti-association df. Defaults to None.
        drop (string): Whether to drop values 'higher' or 'lower' than the threshold. Defaults to 'higher'.
        keep (string): In symmetric matrices, comparing A-B to B-A: keep 'lower', 'higher', or 'both'? Defaults to 'lower'.
        debug (bool, optional): True or string 'deep' for deep debugging. Defaults to False.

    Returns:
        pandas.DataFrame: An edgelist with source, target, value, and anti flag
    """

    if drop not in ['higher', 'lower']:
        raise Exception(f"ERROR! Drop policy should either be 'higher' or lower'. Received value: {drop}")

    if keep not in ['higher', 'lower', 'both']:
        raise Exception(f"ERROR! Keep policy for symmetric matrices should either be 'higher', lower', or 'both'. Received value: {keep}")

    name = None
    if 'name' in data.keys():
        if debug: print(f"Filtering {data['name']}...")
        name = data['name']
    else:
        warnings.warn("This dataframe is unnamed. Please name it to avoid confusion.")
        if debug: print("Filtering unnamed dataframe...")

    df = data['dataframe']

    if selected is None:
        if debug: print(f"No focus column selected. Not filtering {name}.")
    else:
        if debug: print(f"selected: {selected}")
        if selected not in df.columns:
            if selected in df.index:
                # Test to make sure the df isn't accidentally pivoted; fix if so
                # Why? Focus should end up in the middle column of a tripartite network
                df = df.T
            else:
                if debug: print(df.index)
                if debug: print(df.columns)
                raise Exception(f"ERROR! selection {selected} not found in index or columns for {name}. Is this dataset relevant?")

    # Get the total number of values in the df (this will be df.shape[0] x df.shape[1] - nulls)
    original_size = df.count().sum()

    # This if/else only applies if value param is None
    if (value is None) and ('threshold' in data.keys()):
        # Check if threshold was sent in metadata
        value = data['threshold']
        if debug: print(f"Threshold: {value}")
    elif value is None:
        # Set to maximum value
        value = max(df.max())
        if debug: print(f"Threshold: {value} (max, no thresholding applied)")

    # This if/else only applies if anti param is None
    if (anti is None) and ('anti' in data.keys()):
        # Check if anti flag was send in metadata
        anti = data['anti']
        if debug: print(f"Anti: {anti}")
    elif anti is None:
        # Set to default value (assume association)
        anti = False
        if debug: print(f"Anti: {anti} (default)")

    if selected is not None:
        if debug: print(f"Filtering {name} matrix by {selected} {drop} than {value}...")
        lowest = df[selected].min()
        lowest_name = df[selected].idxmin()

        # Drop irrelevant row
        df = df.dropna(subset=[selected])

        # Filter by value
        if drop == 'higher':
            df = df[df[selected] < value]
        elif drop == 'lower':
            df = df[df[selected] > value]
        if df.shape[0] == 0:
            warnings.warn(f"ERROR! No edges found for {selected} with value < {value} in {name}.\nThe lowest value I see is {lowest} for {lowest_name}.")

        # Drop irrelevant cols
        df = df.dropna(axis=1, how='all')

        if debug: print(f"Original edge count: {original_size}\nReduced to {df.count().sum()} edges")

    # Remove identical flipped edges:
    # Check all columns:
    for col in df.columns:
        # If that column is also a row:
        if col in df.index:
            # For each row:
            for row in df.index:
                # If A-B > B-A and keep lower, or A-B < B-A and keep higher, drop A-B
                if ((df.loc[row, col] >= df.loc[col, row]) and (keep == 'lower'))\
                    or ((df.loc[row, col] <= df.loc[col, row]) and (keep == 'higher')):
                        df.loc[row, col] = None
                # If A-B < B-A and keep lower, or A-B > B-A and keep higher, drop B-A
                if ((df.loc[row, col] <= df.loc[col, row]) and (keep == 'lower'))\
                    or ((df.loc[row, col] >= df.loc[col, row]) and (keep == 'higher')):
                        df.loc[col, row] = None
                # Otherwise, keeps both A-B and B-A

    # Convert the matrix to an edge list
    edgelist = df.stack().reset_index().dropna(how='any')
    edgelist.columns = ['source', 'target', 'value']
    if debug: print(f"edgelist shape: {edgelist.shape}\n{edgelist.head()}")

    # Create a boolean flag for anti-associated edges
    edgelist['anti'] = anti

    # Save the source
    edgelist['provenance'] = name

    # Filter edges by value
    # Why filter twice? Stack runs faster with fewer cols
    # but some cols are only partially empty
    # so have to be thresheld after stack
    if drop == 'higher':
        edgelist = edgelist[edgelist['value'] < value]
    elif drop == 'lower':
        edgelist = edgelist[edgelist['value'] > value]

    if selected is not None:
        # Add direct flag
        edgelist['direct'] = [True if selected in [source, target] else False for source, target in zip(edgelist['source'], edgelist['target'])]
    else:
        # When not filtering on a selected column, all edges are 'direct' (will be shown as solid lines)
        edgelist['direct'] = True

    # Calculate the minimum and maximum values of 'value'
    value_min = edgelist['value'].min()
    value_max = edgelist['value'].max()

    # Map the 'value' values to the range of 0.1 to 1
    edgelist['weight'] = edgelist['value'].apply(lambda x: 0.5 + 0.5 * ((x - value_min) / (value_max - value_min)))

    # Cleanup
    edgelist = edgelist.reset_index(drop=True)

    # Store metadata
    edgelist_datalist = {'name': name,
                     'selected': selected,
                     'value': value,
                     'anti': anti}

    edgelist_datalist['dataframe']=edgelist
    return edgelist_datalist

def get_edgelists(dataframe_datalist, selected=None, value=None, debug=False):
    """Take an array of dataframes and return an array of edgelists (basically, run filter_df on all of the dataframes)

    Args:
        dataframes (array): a list of Pandas dataframes to convert
        selected (string, optional): The element of interest to use for filtering. Defaults to None.
        value (float, optional): Override threshold (usually, this should be None and use the df.threshold attribute). Defaults to None.
        debug (bool, optional): Print debugging statements. Defaults to False.

    Returns:
        edgelists: An array of Pandas dataframes of edges with ['source', 'target', 'value', 'weight', 'anti']
    """
    edgelist_datalist = []
    for df_data in dataframe_datalist:
        print(df_data['name'])
        edgelist_data = filter_df(df_data, selected=selected, value=value, debug=debug)
        # Test
        nan_count = edgelist_data['dataframe'].isna().sum().sum()
        if nan_count != 0:
            warnings.warn(f"Heads-up: There are {nan_count} NaNs in the edgelist for {df_data['name']} (Non-fatal, possibly expected.)")
        edgelist_datalist.append(edgelist_data)

    return edgelist_datalist


# %%
def get_nodes(edgelist_datalist, selected=None, debug=False, width=None, height=None, join_substrings=['ENSMUSG', 'g__', 'ranknorm']):
    """Determine the node locations for a graph

    Args:
        edgelist_datalist (array): An array of dicts including key, val 'dataframe': Pandas dataframe of columns source, target, value, weight, anti
        debug (bool or string, optional):  Display print statements if true. If 'deep', display extra print statements for deep debugging. Defaults to False.
        width (int, optional): The width of the chart. Defaults to 1000.
        height (int, optional): The height of the chart. Defaults to 1000.
        join_substrings (list, optional): A list of substrings to join on. Defaults to ['ENSMUSG', 'g__', 'ranknorm'].

    Returns:
        pandas.DataFrame: A dataframe of nodes with locations ['x'] and ['y']
    """

    if width is None: width = 1000
    if height is None: height = 1000

    # This portion checks node sets and merges them if the sets overlap
    # For example, if there are node sets [F, G, H] and [G, H, I],
    # they will be merged into [F, G, H, I]

    # Create the node_sets dict
    node_sets = {}
    set_id=0
    # Iterate over each dataframe in edgelists
    for edgelist_data in edgelist_datalist:
        edgelist = edgelist_data['dataframe']
        # Check 'source' and 'target' columns separately against 'node_sets'
        for column in ['source', 'target']:
            column_values = edgelist[column].tolist()
            if debug == 'deep': print(f"Checking {edgelist.name}: {column} with {len(column_values)} values.")
            new_set = set()
            seen = False

            if len(node_sets)>0:
                if debug: print("There are already sets in the node_sets dict.")
                # Check if any item in the column is already in an existing set
                for key, node_set in node_sets.items():
                    if any(value in list(node_set) for value in column_values):
                        if debug == 'deep': print(f"Found a match in set_id {key} which has a length of {len(node_set)}")
                        # Pull the existing set
                        new_set = node_sets[key]
                        # Add the whole column to the existing set
                        new_set.update(column_values)
                        if debug: print(f"new_set length: {len(new_set)}")
                        # Push back to the node_sets dict
                        node_sets[key] = new_set
                        if debug: print(f"Number of node sets: {len(node_sets)}")
                        seen = True
            # If we need to create a new set
            if seen==False:
                if debug: print(f"Creating a new set at set_id {set_id}")
                # Create a set from the column and push it to the node_sets dict
                node_sets[set_id] = set(column_values)
                if debug == 'deep': print(f"Added set {node_sets[set_id]}")
                if debug: print(f"Number of node sets: {len(node_sets)}")
                # Increment set_id
                set_id+=1
                if debug: print(f"Incremented set_id: {set_id}")


    # Print the resulting node_sets
    for set_id, node_set in node_sets.items():
        if debug == 'deep': print(f"set_id {set_id}: {node_set}")

    # This part merges node sets that share a certain substring
    # For example, if you have node sets [Fx, Gy, Hz] and [Aw, Bx, Cy],
    # and add substring 'x' to join_substrings, sets will be merged to [Fx, Gy, Hz, Aw, Bx, Cy]
    # but adding substring 'z' to join_substrings will not merge the sets

    # Check all the specified substrings for groupings
    for substring in join_substrings:
        if debug == 'deep': print(f"Looking for substring '{substring}' in multiple sets")
        to_join = []
        to_remove = []
        for set_id, node_set in node_sets.items():
            for n in node_set:
                if substring in n:
                    # Save to join list the first time you find a matching element
                    to_join.append(node_set)
                    to_remove.append(set_id)
                    if debug == 'deep': print(f"Found '{substring}' in set {set_id} in {n}")
                    break
        # If we found things to join, join them
        if len(to_join)>1:
            if debug == 'deep': print(f"Found {len(to_join)} sets to join")
            new_set = set.union(*to_join)
            for i in to_remove:
                # Take out the subsets
                node_sets.pop(i)
                if debug == 'deep': print(f"Popped set {i}")
            # Add the new superset
            set_id+=1
            node_sets[set_id] = new_set
            if debug == 'deep': print(f"Created a new set with set_id {set_id}")
            set_id +=1
            for key, node_set in node_sets.items():
                if debug == 'deep': print(f"set_id {key}: {node_set}")
            if debug == 'deep': print("\n")

    # Delete key: set pairs with length 0
    node_sets = {key: values for key, values in node_sets.items() if len(values) > 0}

    # This part assigns layers to the nodes, and places the selected node in the middle

    # Create an empty DataFrame
    nodelist = pd.DataFrame(columns=['layer'])
    layer_count = len(node_sets.keys())
    i=0
    middle_layer = layer_count // 2
    # Iterate over the dictionary items
    for layer, values in node_sets.items():
        if len(values)>0:
            values = list(values)
            if debug: print(f"Placing layer that begins with {values[0]}. i={i}.")
            old_i = i
            is_selected = False
            # Check if the selected node is in this layer
            # If it is, we'll need to place this layer horizontally in the middle
            # and move the selected node vertically to the top
            # Note: Doesn't break if selected is None
            if selected in values:
                # Position the selected node at the top
                values.remove(selected)
                values.insert(0, selected)
                values.reverse()
                # Set to the index of the middle layer
                i = middle_layer
                if debug: print(f"Found selected item {selected} in {layer}. Temporarily setting i={i} (middle layer).")
                # We need this because we don't increment i if we're on the middle layer
                is_selected = True
            else:
                # If we're on the middle layer but it doesn't contain the selected item,
                # skip ahead one layer before adding
                if i == middle_layer: i += 1
            # Create a temporary DataFrame for each layer
            layer_df = pd.DataFrame({
                'layer': i,
                'node': values,
                'x': (width / layer_count) * (i) # Distribute the layers evenly across the width
                })
            layer_df.reset_index(drop=True, inplace=True)
            # Find the vertical spacer for this layer
            if len(values) > 1:
                # If there's more than one node, divide the height by the number of nodes
                spacer = height / (len(values)-1)
            else:
                # If there's only one node in the layer, place it in the vertical middle
                spacer = height / 2
            # Set the y-location
            layer_df['y'] = [(spacer * (j)) for j in layer_df.index]
            # Concatenate the temporary DataFrame to the main DataFrame
            nodelist = pd.concat([nodelist, layer_df])
            # Resume counting; don't add 1 if we're on the middle layer
            i = old_i + 1 if is_selected == False else old_i

    # Reset the index of the final DataFrame
    nodelist.reset_index(drop=True, inplace=True)

    return nodelist

def get_edge_positions(edgelist_datalist, nodelist, debug=False):
    """Merge the nodelist dataframe with the edgelist dataframe and perform necessary transformations.

    Args:
        edgelist (pandas.DataFrame): The edgelist dataframe with columns 'source', 'target', 'value', 'anti'.
        nodelist (pandas.DataFrame): The nodelist dataframe with columns 'layer', 'node', 'x', 'y'.
        debug (bool, optional): Print debugging statements. Defaults to False.

    Returns:
        pandas.pandas.DataFrame: The updated edgelist dataframe with added columns 'x1', 'y1', 'x2', 'y2', 'color', 'dash'.
    """

    # Merge all the edgelists into one dataframe
    edgelists = []
    for edgelist_data in edgelist_datalist:
        edgelists.append(edgelist_data['dataframe'])

    edgelist = pd.concat(edgelists).reset_index(drop=True)
    if edgelist.isna().sum().sum() != 0:
            warnings.warn(f"There are NaNs in the edgelist when combining in get_edge_positions!")

    # Merge the nodelist dataframe with the edgelist dataframe for source nodes
    edgelist = edgelist.merge(nodelist, left_on='source', right_on='node', how='left')

    # Rename the x and y columns for the source nodes
    edgelist = edgelist.rename(columns={'x': 'x1', 'y': 'y1'})

    # Drop the unnecessary columns from the source merge
    edgelist = edgelist.drop(['node', 'layer'], axis=1)

    # Merge the nodelist dataframe with the edgelist dataframe for target nodes
    edgelist = edgelist.merge(nodelist, left_on='target', right_on='node', how='left')

    # Rename the x and y columns for the target nodes
    edgelist = edgelist.rename(columns={'x': 'x2', 'y': 'y2'})

    # Drop the unnecessary columns from the target merge
    edgelist = edgelist.drop(['node', 'layer'], axis=1)

    return edgelist

def format_text(text):
    """Format text for display in the network graph. Replace underscores with spaces
    and split long strings into two lines.

    Args:
        text (string): The text to format.

    Returns:
        string: The formatted text.
    """
    if '_' in text:
        text = ' '.join(text.split('_')[1:])
    if (len(text)>20) and (' ' in text):
        matches = re.finditer(' ', text)
        match_indices = [match.start() for match in matches]
        split_spot = match_indices[len(match_indices) // 2]
        text = text[:split_spot] + '<br>' + text[split_spot:]
    return text
# %%
