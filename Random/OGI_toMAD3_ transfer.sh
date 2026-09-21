#!/bin/bash

tmux attach -t mad3_transfer2
tmux new -s mad3_transfer2

# Inside the tmux session:

Destination=/gpfs/fs1/home/jl2251/mount/data/users/
source_path=/data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz

/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > DoubleRsync_MAD3_transfer_yluoTar1.log 2>&1 &

source_path=/data/Segre_Lab/users/jlama/yluo_archive_7July2026.log

/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > rsync_MAD3_transfer_yluoTar2.log 2>&1 &

wait
echo "Both transfers finished"

you can just detach (Ctrl+b, then d)

tmux attach -t mad3_transfer2
jobs                          # shows which background jobs are still running
tail -20 rsync_MAD3_transfer_yluoTar1.log
tail -20 rsync_MAD3_transfer_yluoTar2.log

######### SPlit the user directory tar ball in half and then transfer
tmux attach -t mad3_transfer2
source_path=/data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz
split -n 2 ${source_path} ${source_path}.part_ > splittedTar.log 2>&1 &
#### Job running as of Aug 27 2026 #### 
#rm "${source_path}"

Destination=/gpfs/fs1/home/jl2251/mount/data/DRCR/

/usr/bin/rsync -ravoP "${source_path}.part_aa" "${Destination}" > rsync_MAD3_transfer_yluoTar_partaa.log 2>&1 &
/usr/bin/rsync -ravoP "${source_path}.part_ab" "${Destination}" > rsync_MAD3_transfer_yluoTar_partab.log 2>&1 &

#Double RSYNC
#/usr/bin/rsync -ravoP "${source_path}.part_aa" "${Destination}" > Doublersync_MAD3_transfer_yluoTar_partaa.log 2>&1 &
#/usr/bin/rsync -ravoP "${source_path}.part_ab" "${Destination}" > Doublersync_MAD3_transfer_yluoTar_partab.log 2>&1 &


#### Job running as of Sept 21 2026 ####  DONE

cd ${Destination}
cat users_yluo_7July2026.tar.gz.part_aa users_yluo_7July2026.tar.gz.part_ab > users_yluo_7July2026.tar.gz 2> tarball.merge.log &

#Original main tar file
ef3d3af37b3262ed7c1ef9dce45ee09d7710cc0763c2db57e115c1f3c3b0b778  /data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz

#### Job rerunning as of Aug 31 2026 #### Done
#Original split tar file
sha256sum /data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz.part_aa > /data/Segre_Lab/users/jlama/users_yluo_7July2026_aa.tar.gz.OGI.sha256 2> /data/Segre_Lab/users/jlama/sha256.yluo_aa_error.log &
sha256sum /data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz.part_ab > /data/Segre_Lab/users/jlama/users_yluo_7July2026_ab.tar.gz.OGI.sha256 2> /data/Segre_Lab/users/jlama/sha256.yluo_ab_error.log &

e0dcae71a21b1e48bd1c40a3324d6ee50bff5186e845edbc4043b70fcae3f763  /data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz.part_aa
58d8cc520df3065cb1339a8bf6647ca743a6fa60695c284c0de586ab3bd8ea95  /data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz.part_ab

#Mad3 transferred split tar file
sha256sum users_yluo_7July2026.tar.gz.part_aa > users_yluo_7July2026_aa_des.tar.gz.sha256 2> sha256.yluo_aa_error_dest.log &
sha256sum users_yluo_7July2026.tar.gz.part_ab > users_yluo_7July2026_ab_des.tar.gz.sha256 2> sha256.yluo_ab_error_dest.log &
sha256sum users_yluo_7July2026.tar.gz > users_yluo_7July2026.des.tar.gz.sha256 2> sha256.yluo_errors_2.des.log &

e0dcae71a21b1e48bd1c40a3324d6ee50bff5186e845edbc4043b70fcae3f763  users_yluo_7July2026.tar.gz.part_aa
58d8cc520df3065cb1339a8bf6647ca743a6fa60695c284c0de586ab3bd8ea95  users_yluo_7July2026.tar.gz.part_ab
99aac7b029e738f2a7fd7ebdefa321539dc4bc4143037b6ee9090728c7c46a70  users_yluo_7July2026.tar.gz


diff <(awk '{print $1}' users_yluo_7July2026.des.tar.gz.sha256) \
     <(awk '{print $1}' users_yluo_7July2026.tar.gz.sha256)
     
# Test gzip integrity (checks compressed stream isn't corrupted)
gzip -t users_yluo_7July2026.tar.gz && echo "gzip OK" &
# List tar contents without extracting, to confirm it reads through fully
tar -tzf users_yluo_7July2026.tar.gz > /dev/null && echo "tar structure OK" &
     
#### Job rerunning as of Aug 31 2026 #### Running

#Verify
# Compare size to the original
stat --format="%s" users_yluo_7July2026.tar2.gz
# should match: 4800047099227

# Confirm the reassembled archive is actually intact
gzip -t users_yluo_7July2026.tar2.gz && echo "OK — archive is valid" &

# Only once verified, clean up the two halves
rm users_yluo_7July2026.tar.gz.part_aa users_yluo_7July2026.tar.gz.part_ab users_yluo_7July2026.tar.gz

############################################

tmux new -s mad3_transfer_DRCR
source_path=/gpfs/fs1/data/Segre_Lab/data/DRCR_13July2026.tar.gz
Destination=/gpfs/fs1/home/jl2251/mount/data/DRCR/

/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > rsync_MAD3_transfer_DRCR.log 2>&1 &
/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > rsync_MAD3_transfer_DRCR2.log 2>&1 &

wait
echo "Both transfers finished"

you can just detach (Ctrl+b, then d)

tmux attach -t mad3_transfer_DRCR
jobs                          # shows which background jobs are still running
tail -20 rsync_MAD3_transfer_DRCR.log
tail -20 rsync_MAD3_transfer_DRCR2.log
#############################################


echo starting transfer of SCORE tar file
source_path=/gpfs/fs1/data/Segre_Lab/data/SCORE_13July2026.tar.gz
#/usr/bin/rsync -ratlzv --rsh="$HOME/.conda/envs/mytools/bin/sshpass -p $password ssh -o StrictHostKeyChecking=no -l jl2251" eris2n8.research.partners.org:${source_path}  ${Destination} ; echo "Exit code: $?" >> transfer.log

source_path=/gpfs/fs1/data/Segre_Lab/data/SCORE_archive_13July2026.log
#/usr/bin/rsync -ratlzv --rsh="$HOME/.conda/envs/mytools/bin/sshpass -p $password ssh -o StrictHostKeyChecking=no -l jl2251" eris2n8.research.partners.org:${source_path}  ${Destination} ; echo "Exit code: $?" >> transfer.log

echo Done!
  
  #Round 1
  
  ######################
#!/bin/bash

Destination=/gpfs/fs1/home/jl2251/mount/data/users/
source_path=/data/Segre_Lab/users/jlama/users_yluo_7July2026.tar.gz

for i in 1 2 3 4 5; do
  rsync -ravoP --partial --partial-dir=.rsync-partial \
    --timeout=300 \
    "${source_path}" "${Destination}" \
    >> rsync_MAD3_transfer_yluoTar1.log 2>&1
  status=$?
  [[ $status -eq 0 ]] && { echo "Success on attempt $i"; break; }
  echo "Attempt $i failed (exit $status), retrying in 60s..."
  sleep 60
done


#############
tmux attach -t mad3_transfer_Sobrin_R01
tmux new -s mad3_transfer_Sobrin_R01
./Sobrin_R01_tar.sh > sobrinR01.tar.log 2>&1 &

Destination=/gpfs/fs1/home/jl2251/mount/data/
source_path=/gpfs/fs1/data/Segre_Lab/data/Sobrin_R01_GC_OHTN.tar.gz.tar.gz

source_path=/gpfs/fs1/data/Segre_Lab/data/Sobrin_R01_GC_OHTN.tar.gz_27Aug2026.log
/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > rsync_MAD3_transfer_Sobrin_R01_GC_OHTN_2.log 2>&1 &

/usr/bin/rsync -ravoP "${source_path}" "${Destination}" > Rersync_MAD3_transfer_Sobrin_R01_GC_OHTN_1.log 2>&1 &
sha256sum "${source_path}" > /gpfs/fs1/data/Segre_Lab/data/Sobrin_R01_GC_OHTN.tar.gz.sha256 2> sha256.Sobrin_R01_GC_OHTN_errors.log &
sha256sum ${Destination}/Sobrin_R01_GC_OHTN.tar.gz.tar.gz > ${Destination}/Sobrin_R01_GC_OHTN.tar.gz.sha256 2> sha256.Sobrin_R01_GC_OHTN_errors.destin.log &

#Original main tar file
55fcd3efd1213d52043e0f93b9da6071ccb29f02e6a1ef05b1a338d08c81ab4c  /gpfs/fs1/data/Segre_Lab/data/Sobrin_R01_GC_OHTN.tar.gz.tar.gz
55fcd3efd1213d52043e0f93b9da6071ccb29f02e6a1ef05b1a338d08c81ab4c  /gpfs/fs1/home/jl2251/mount/data//Sobrin_R01_GC_OHTN.tar.gz.tar.gz

#TRUE
rm -r /gpfs/fs1/data/Segre_Lab/data/Sobrin_R01_GC_OHTN.tar.gz.tar.gz
#### Job running as of Aug 31 2026 #### Done 


###################### Confirm Tar files transfer ################################


SOURCE=/gpfs/fs1/data/Segre_Lab/data/
tar -czvf ${SOURCE}/SCORE_13July2026.tar.gz -C ${SOURCE} SCORE > ${SOURCE}/SCORE_archive_13July2026.log 

SOURCE=/gpfs/fs1/data/Segre_Lab/data/
tar -czvf ${SOURCE}/DRCR_13July2026.tar.gz -C ${SOURCE} DRCR > ${SOURCE}/DRCR_archive_13July2026.log 

########### For DRCR
source_path=/gpfs/fs1/data/Segre_Lab/data/DRCR_13July2026.tar.gz
Destination=/gpfs/fs1/home/jl2251/mount/data/DRCR/

# BEFORE moving, on the source server:
sha256sum "${source_path}" > /gpfs/fs1/data/Segre_Lab/data/DRCR_13July2026.tar.gz.sha256 2> sha256_errors.log &

sha256_pid=$!
echo "sha256sum running as PID ${sha256_pid}"
if kill -0 "${sha256_pid}" 2>/dev/null; then
  echo "Still running..."
else
  echo "Finished."
fi
# AFTER moving, on the destination server:
cd ${Destination}
sha256sum -c ${Destination}/DRCR_13July2026.tar.gz > DRCR_13July2026.tar.gz.sha256 2>&1 &

########### For SCORE
source_path=/gpfs/fs1/data/Segre_Lab/data/SCORE_13July2026.tar.gz
Destination=/gpfs/fs1/home/jl2251/mount/data/SCORE/

# BEFORE moving, on the source server:
sha256sum "${source_path}" > /gpfs/fs1/data/Segre_Lab/data/SCORE_13July2026.tar.gz.sha256 2> sha256.score_errors.log &

sha256_pid=$!
echo "sha256sum running as PID ${sha256_pid}"
if kill -0 "${sha256_pid}" 2>/dev/null; then
  echo "Still running..."
else
  echo "Finished."
fi
# AFTER moving, on the destination server:
cd ${Destination}
sha256sum -c ${Destination}/SCORE_13July2026.tar.gz > SCORE_13July2026.tar.gz.sha256 2>&1 &

############################################################################################


### Checking validity of transferred tar file ####

# 1. Check exact file size matches source
ls -la /gpfs/fs1/home/jl2251/mount/data/users/users_yluo_7July2026.tar.gz

# Compare against source size: 4,800,047,099,227 bytes

md5sum users_yluo_7July2026.tar.gz > sourceFile_users_yluo.md5sum.txt 2>&1 &
md5sum /gpfs/fs1/home/jl2251/mount/data/users/users_yluo_7July2026.tar.gz > destinationFile_users_yluo.md5sum.txt 2>&1 &

# 2. Test archive integrity — confirms it's not truncated or corrupted
gzip -t /gpfs/fs1/home/jl2251/mount/data/users/users_yluo_7July2026.tar.gz && echo "OK: archive intact"

# 3. List contents of the tar without extracting — will error out if truncated
tar -tzf /gpfs/fs1/home/jl2251/mount/data/users/users_yluo_7July2026.tar.gz > /dev/null && echo "OK: tar structure intact"

#### Sept 21 2026 #####

SOURCE=/gpfs/fs1/data/Segre_Lab/data/
tar -czvf ${SOURCE}/DRCR_21Sept2026.tar.gz -C ${SOURCE} DRCR > ${SOURCE}/DRCR_archive_21Sept2026.log 2> DRCR_tar.log &
sha256sum "${SOURCE}/DRCR_21Sept2026.tar.gz" > "${SOURCE}/DRCR_21Sept2026.tar.gz.sha256" 2> "${SOURCE}/DRCR_21Sept2026_error.log" &

### Current running jobs
tmux attach -t mad3_transfer2
[1]-  Running                 cat users_yluo_7July2026.tar.gz.part_aa users_yluo_7July2026.tar.gz.part_ab > users_yluo_7July2026.tar.gz 2> tarball.merge.log &  (wd: ~/mount/data/users)
[2]+  Running                 tar -czvf ${SOURCE}/DRCR_21Sept2026.tar.gz -C ${SOURCE} DRCR > ${SOURCE}/DRCR_archive_21Sept2026.log 2> DRCR_tar.log &
