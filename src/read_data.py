import numpy as np


def read_data_solution(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    u_00 u_01 ... u_0n
    u_10 u_11 ... u_1n
    ...
    u_n0 u_n1 ... u_nn

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data


def read_data_file(file_path):
    headers = []
    data_blocks = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line:
            header = float(line)
            headers.append(header)

            i += 1
            block_lines = []
            while i < len(lines) and lines[i].strip():
                aux = list(map(float, lines[i].split()))
                # repeat the first element of each line to make the
                # animation periodic
                aux.append(aux[0])
                block_lines.append(aux)
                i += 1
            # repeat the first line to make the animation periodic
            block_lines.append(block_lines[0])
            data_block = np.array(block_lines)
            data_blocks.append(data_block)
        i += 1

    headers_array = np.array(headers)
    data_blocks_array = np.array(data_blocks)

    # now we dupplicate the

    return headers_array, data_blocks_array


def read_data_file_freq(file_path):
    headers = []
    data_blocks = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line:
            header = float(line)
            headers.append(header)

            i += 1
            block_lines = []
            while i < len(lines) and lines[i].strip():
                aux = list(map(float, lines[i].split()))
                block_lines.append(aux)
                i += 1
            data_block = np.array(block_lines)
            data_blocks.append(data_block)
        i += 1

    headers_array = np.array(headers)
    data_blocks_array = np.array(data_blocks)

    # now we dupplicate the

    return headers_array, data_blocks_array


def read_data_energy(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    t_1 E_1 dE_1 u(pi, pi) x_max1 y_max1
    t_2 E_2 dE_2 u(pi, pi) x_max2 y_max2
    ...

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data


def read_data_energy_return(filename: str) -> np.ndarray:
    """Read data from filename. The data is of the form:
    t1 E_1
    t2 E_2
    ...

    Args:
        filename (str): File path.

    Returns:
        np.ndarray: Data.
    """
    data = np.loadtxt(filename)
    return data
