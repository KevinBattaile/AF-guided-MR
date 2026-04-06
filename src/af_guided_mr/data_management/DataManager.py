import gemmi
import logging

class DataManager:
    def __init__(self) -> None:
        pass

    def get_space_group(self, mtz_path):
        mtz = gemmi.read_mtz_file(mtz_path)
        space_group = mtz.spacegroup
        return space_group.hm

    def get_high_resolution(self, mtz_path):
        mtz = gemmi.read_mtz_file(mtz_path)
        return mtz.resolution_high()

    @staticmethod
    def get_mtz_labels(mtz_path):
        """
        Scans MTZ file to preemptively select the best amplitude/intensity and Free R flag columns.
        Excludes SA_flag to prevent ambiguous array issues in Phenix.
        """
        mtz = gemmi.read_mtz_file(mtz_path)
        columns = [col.label for col in mtz.columns]

        # Priorities for Data labels
        priorities = [
            ("IMEAN", "SIGIMEAN"),
            ("I", "SIGI"),
            ("F", "SIGF"),
            ("FP", "SIGFP"),
        ]

        selected_data_labels = None
        for prio in priorities:
            if prio[0] in columns and prio[1] in columns:
                    selected_data_labels = f"{prio[0]},{prio[1]}"
                    break

        # Priorities for Free R labels
        free_r_priorities = ["FreeR_flag", "R-free-flags", "FREE"]
        selected_free_r_label = None
        for prio in free_r_priorities:
            if prio in columns and prio != "SA_flag":
                selected_free_r_label = prio
                break

        if selected_data_labels is None:
            logging.warning("Could not find suitable data labels in MTZ file. Refinement may fail.")
        if selected_free_r_label is None:
            logging.warning("Could not find suitable Free R label in MTZ file. Refinement may fail.")

        return selected_data_labels, selected_free_r_label
