git clone https://github.com/deepomicslab/SpecHLA.git --depth 1
cd SpecHLA/
conda env create --prefix=./spechla_env -f environment.yml
conda activate ./spechla_env
# Installing a Printer in Linux

Here's how to install a printer on Linux systems:

## Using CUPS Web Interface (Works on all distributions)

1. Make sure CUPS is installed:
   ```bash
   sudo apt install cups   # For Debian/Ubuntu
   sudo dnf install cups   # For Fedora/RHEL
   ```

2. Start/enable the CUPS service:
   ```bash
   sudo systemctl start cups
   sudo systemctl enable cups
   ```

3. Access the CUPS web interface:
   - Open a browser and navigate to: `http://localhost:631`
   - Go to "Administration" > "Add Printer"
   - Enter your username and password when prompted
   - Follow the wizard to select your printer and install drivers

## Using Desktop Environment Tools

### GNOME (Ubuntu, Fedora):
```bash
sudo apt install system-config-printer   # On Debian/Ubuntu
```
1. Open "Settings" > "Printers"
2. Click "Add Printer"
3. Select your printer from the list or use "Find Network Printer"

### KDE:
1. Open "System Settings" > "Printers"
2. Click "Add" > "Add Printer"

## Command Line Method

```bash
# List available printers
lpinfo -v

# Add a network printer (example)
lpadmin -p PrinterName -E -v ipp://printer.local/ipp/print -m everywhere
```

## Installing Printer Drivers

For many printers, you may need specific drivers:
```bash
sudo apt install printer-driver-all    # General drivers package
sudo apt install hplip                 # For HP printers
```

After installation, print a test page to verify everything works correctly.