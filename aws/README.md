# Running Optimization Jobs on AWS with Starcluster
Requirement : AWS and Wasabi accounts

#### Create an AMI with necessary software:
* Launch a t2 micro instance and install all the necessary software on a conda environment.
* Set up wasabi: 
```sh
$ pip install awscli-plugin-endpoint
$ aws configure # use access key and secret access keys, default output: json, region us-east-1
$ aws configure --profile wasabi # now use wasabi access keys and secret access keys,region us-west-1
$ aws configure set plugins.endpoint awscli_plugin_endpoint
$ aws configure --profile wasabi set s3.endpoint_url https://s3.us-west-1.wasabisys.com
$ aws s3 ls --profile wasabi  # check wasabi buckets
$ aws s3 cp dir s3://aibs.snmo.01/dir_group/dir --recursive --profile wasabi
```
* Log out

#### Create an EBS volume

* Create an EBS volume using the AWS management console 
* Attach the EBS volume to the running t2 micro instance.
* Log in to the t2 micro instance and mount the EBS volume
```sh
$ lsblk # Generally something like xvdf would show up
$ sudo file -s /dev/xvdf # Get information about filesystem
$ sudo mkfs -t xfs /dev/xvdf # Warning : only for empty volume
$ sudo mkdir /data
$ sudo mount /dev/xvdf /data
```
* Unmout the EBS volume
```sh
$ sudo umount /data
```
* Detach the EBS volume from the instance using AWS management console

#### Create a snapshot of the AMI
* AWS Management Console -> Instances -> Actions -> Image -> Create Image
* Give it a name and increase the volume a little bit while saving the snapshot.


