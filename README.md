# CCorticosteroids for infectious critical illness: A multicenter target trial emulation stratified by predicted organ dysfunction trajectory

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

```
# Clone the repository
git clone https://github.com/surajraj99/Corticosteroids-in-Patients-with-Sepsis.git

# Navigate to the project directory
cd Corticosteroids-in-Patients-with-Sepsis

# Install dependencies (if any)
pip install -r requirements.txt
```

## Usage

### Data Creation and Sample Data

Scripts for creating the data and some sample data for the eICU and MIMIC cohorts are located in the `data` folder. Not all filepaths are correct and need to be rewritten to fit the user's file structure. For subtyping, refer to https://github.com/xuzhenxing2019/sepsis_subphenotype.

### Trial Emulation Scripts

The main scripts for running the trial emulation are located in the main directory. To run the trial emulation for each of the different outcomes, run the releveant notebooks (outcome is in the title). Sample data for running the files are already in the `data` folder. Only data and trial emulation script for eICU-MIMIC are available.

## Directory Structure

```
.
├── 28 Day Mortality.ipynb
├── Time to Discharge.ipynb
├── Cessation of Mechanical Ventilation.ipynb
├── Predict Zhenxing Subtypes - MIMIC.ipynb
├── Predict Zhenxing Subtypes - eICU.ipynb
├── data/
│   ├── eICU
|   |   ├── sample data files
|   |   ├── scripts to create some eICU data
|   |        ├── scripts
│   ├── MIMIC
|       ├── sample data files
|       ├── scripts to create some MIMIC data
|            ├── scripts
│   └── MIMIC
└── README.md
```

## Contributing

If you would like to contribute to this project, please fork the repository and submit a pull request.

## Contact

For more information, feel free to reach out to:

**Suraj Rajendran**  
📧 Email: [sur4002@med.cornell.edu](mailto:sur4002@med.cornell.edu)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Cite Us
Rajendran S, Xu Z, Pan W, Zhang C, Schenck E, Wang F. Corticosteroids in Patients with Sepsis: An Observational Research through Target Trial Emulations
