import yaml
import random
from copy import deepcopy
from pathlib import Path

input_file = "part_params.yaml"
output_dir = Path("sampled_yaml")
n_files = 2   # number of files to generate

# lower/upper uncertainty boundary
uncertainty = {
    #"Fe": (0.20, 0.30),
    "O2": (0.30, 0.30),
    "O": (0.3, 0.3),
    "H2O": (0.3, 0.3),
    "H2": (0.3, 0.3),
    "H": (0.3, 0.3),
    "OH": (0.3, 0.3)
}

def sample_ta_truncated(base_ta, lower_frac, upper_frac, max_tries=1000):
    low = base_ta * (1.0 - lower_frac)
    high = base_ta * (1.0 + upper_frac)

    mu = base_ta
    sigma = (high - low) / 4.0   # approx 95% interval

    if sigma <= 0.0:
        return base_ta

    for _ in range(max_tries):
        x = random.gauss(mu, sigma)
        if low <= x <= high:
            return x

    return min(max(mu, low), high)


# Read input YAML once
with open(input_file, "r") as f:
    base_data = yaml.safe_load(f)
    
output_dir.mkdir(exist_ok=True)

# Generate many sampled files
for k in range(1, n_files + 1):
    new_data = deepcopy(base_data)

    if "collisions" in new_data:
        for coll in new_data["collisions"]:
            species = coll.get("name", None)

            if species in uncertainty and "T_a" in coll:
                base_ta = float(coll["T_a"])
                lower_frac, upper_frac = uncertainty[species]

                new_ta = sample_ta_truncated(base_ta, lower_frac, upper_frac)
                coll["T_a"] = float(new_ta)

    out_file = output_dir / f"collisions_sampled_{k:03d}.yaml"

    with open(out_file, "w") as f:
        yaml.dump(new_data, f, sort_keys=False)

print(f"Generated {n_files} files in: {output_dir}")
