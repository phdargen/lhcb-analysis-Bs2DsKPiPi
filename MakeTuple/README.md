# Collect list of LFNs (output stored in LFNs.txt):
# (on lxplus)

ganga -i collect-LFNs.py jobnumber

# Download files (modify outputLocation in collect-output.sh):

ssh -p28 sigma0
source collect-output.sh LFNs.txt



