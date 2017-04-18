## ECEN773 Project Writeup

This directory contains files related to the final project writeup for Jeff's ECEN773 course. The folder pattern of "latex_<name>" can be used for other writeup scenarios as well.

1. Install LaTeX:
    ```bash
    sudo apt install -y texlive-full
    # if you have limited space:
    sudo apt purge --auto-remove texlive-lang-*
    sudo apt purge --auto-remove texlive-*-doc

    ```

1. Compile the .tex source:
    ```bash
    pdflatex writeup.tex
    ```
