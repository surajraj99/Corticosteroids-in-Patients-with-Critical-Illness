# Corticosteroids in Patients with Sepsis: An Observational Research through Target Trial Emulations

## Abstract Summary

This retrospective multicenter study employs a target trial emulation framework to evaluate the effectiveness of corticosteroids in sepsis patients. The study reveals that corticosteroids' effectiveness varies based on machine learning-based stratification, offering nuanced insights for tailored treatment strategies.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Contributing](#contributing)
- [Contact](#contact)
- [License](#license)

## Installation

```bash
# Clone the repository
git clone https://github.com/surajraj99/Corticosteroids-in-Patients-with-Sepsis.git

# Navigate to the project directory
cd Corticosteroids-in-Patients-with-Sepsis

# Install dependencies (if any)
pip install -r requirements.txt
```

## Usage

### Trial Emulation Scripts

The main scripts for running the trial emulation are located in the main directory. To run the trial emulation, execute:

\`\`\`bash
python main_script.py
\`\`\`

### Data Creation and Sample Data

Scripts for creating the data and some sample data for the eICU and MIMIC cohorts are located in the `data` folder.

\`\`\`bash
# Navigate to the data folder
cd data

# Run the data creation script
python create_data.py
\`\`\`

## Directory Structure

\`\`\`
.
├── main_script.py
├── data/
│   ├── create_data.py
│   └── sample_data.csv
└── README.md
\`\`\`

## Contributing

If you would like to contribute to this project, please fork the repository and submit a pull request.

## Contact

For more information, feel free to reach out to:

**Suraj Rajendran**  
📧 Email: [sur4002@med.cornell.edu](mailto:sur4002@med.cornell.edu)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
