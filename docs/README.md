## Troubleshooting

### System
* I receive the error ```error from nextflow (ERR): mkfifo(/tmp/17.inpipe1) failed.``` when using the latest version of nextflow ```nextflow/23.04.3```. Changing to an earlier version removes the error.

### Containers

* In case you receive the error, it is possible that you don't have enough loop devices or they are not specified correctly
```
container creation failed: mount /proc/self/fd/3->/usr/local/var/singularity/mnt/session/rootfs error: while mounting image /proc/self/fd/3: failed to find loop device: could not attach image file to loop device: no loop devices available
```
Some useful commands
* Check availability of loop devices
```
losetup -f
```
If you don't see the loops available, you can try     
```
sudo mknod -m660 /dev/loop8 b 7 8
```
Otherwise you need to specify the number of loops to use by editing ```/etc/modprobe.d/loop.conf```, add the following line to the file      
```
options loop max_loop=100
```
Remember to reload the operation and reinsert it with     
```
sudo rmmod loop
sudo modprobe loop
```
After this process restart the singularity device. I usually do a reboot of the VM.    
