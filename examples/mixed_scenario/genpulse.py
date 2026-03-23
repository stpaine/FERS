import numpy as np
import h5py
import os

"""
To use this script first install numpy and h5py using pip:
    pip install numpy h5py
"""


def generate_pulse_file(filename="pulse.h5", sample_rate=10e6, pulse_width=1e-6):
    """
    Generates an HDF5 file containing a simple rectangular pulse for FERS.
    The file structure matches what libfers' hdf5_handler expects.

    At 10 MHz sample rate and 1 us pulse width, this generates exactly 10 samples.
    """
    num_samples = int(sample_rate * pulse_width)
    i_data = np.ones(num_samples, dtype=np.float64)
    q_data = np.zeros(num_samples, dtype=np.float64)

    if os.path.exists(filename):
        os.remove(filename)

    with h5py.File(filename, 'w') as f:
        # Create I component group and dataset
        i_group = f.create_group("I")
        i_group.create_dataset("value", data=i_data)

        # Create Q component group and dataset
        q_group = f.create_group("Q")
        q_group.create_dataset("value", data=q_data)

    print(f"Generated pulse file '{filename}' with {num_samples} samples.")


if __name__ == "__main__":
    generate_pulse_file()
