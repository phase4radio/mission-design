# Phase4 Link Budget Tool

This folder contains a Python based Link Budget Tool designed to be used interactively with the Jupyter Notebook.

This tool differs from the classic spreadsheet approach in a few ways:

 1. The satellite, channel and groundstations are defined in YAML configuration files for quickly arranging analysis of different scenarios.
 2. The tool dynamically calculates the path loss from the ITU models to get accurate data based on frequency, location and probability.
 3. The tool will dynamically find the best DVB-S2X MODCOM as the link conditions change
 4. The tool supports Monte-Carlo simulations of the link to get statistical insight into the nature of the link

This work merely arranges data around these excellent pieces of work:
 - (pylink)[https://github.com/harrison-caudill/pylink]
 - (ITU-Rpy)[https://github.com/iportillo/ITU-Rpy]
 
