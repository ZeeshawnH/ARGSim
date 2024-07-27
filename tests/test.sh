set -xe
./labp_v22
diff outLABP_original.sites outLABP_[0-9]*.sites
diff outLABP_original.stats outLABP_[0-9]*.stats
rm outLABP_[0-9]*.sites
rm outLABP_[0-9]*.stats