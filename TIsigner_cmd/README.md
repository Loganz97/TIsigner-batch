# TIsigner Command Line Version (Modified for Batch Processing in Google Colab)

This is a modified version of TIsigner, optimized for batch processing of DNA sequences using Google Colab.

### Google Colab Usage

This version is designed to be used with Google Colab. A Colab notebook (`TIsigner_batch.ipynb`) is provided in the repository for easy use.

To use TIsigner in Google Colab:

1. Open the `TIsigner_batch.ipynb` notebook in Google Colab.
2. Run the cells in order. The notebook will:
   - Set up the Conda environment
   - Clone the TIsigner-batch repository
   - Install necessary dependencies

3. Prepare your input CSV file with columns: 'Sequence Name', 'Sequence', and 'Optimized Sequence'.
4. Upload your CSV file when prompted by the notebook.
5. The script will process each sequence and output a consolidated CSV with the results.

### Output

The output is a consolidated CSV file containing:

- Sequence Name
- Original Sequence
- Input Sequence (pre-optimized)
- TISigned Sequence (further optimized by TIsigner)
- Opening Energy
- Score (if available)

### Dependencies

The necessary dependencies are installed automatically when running the Colab notebook. The main requirements are:

- Python 3.6+
- ViennaRNA suite
- Pandas
- NumPy

### Acknowledgments

This work builds upon the original TIsigner project by [Bikash Kumar Bhandari](https://bkb3.github.io), [Chun Shem Lim](https://github.com/lcscs12345), and [Paul P Gardner](https://github.com/ppgardne) (2019-).

### Funding

This work was conducted as part of a <b><font color='green'>US Department of Energy SCGSR Fellowship</font></b> ([details](https://science.osti.gov/wdts/scgsr)) with additional support from the <b><font color='DodgerBlue'>US National Science Foundation</font></b> ([grant](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2132183&HistoricalAwards=false)).
