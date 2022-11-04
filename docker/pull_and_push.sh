#!/usr/bin/env bash

# This script shows how to build the Docker image and push it to ECR to be ready for use
# by SageMaker.

# The argument to this script is the image name. This will be used as the image on the local
# machine and combined with the account and region to form the repository name for ECR.
image=$1
version=$2

if [ "$image" == "" ]
then
    echo "Usage: $0 <image-name> <image-version>"
    exit 1
fi

if [ "$version" == "" ]
then
    echo "Usage: $0 <image-name> <image-version>"
    exit 1
fi


# Get the account number associated with the current IAM credentials
account=$(aws sts get-caller-identity --query Account --output text)

if [ $? -ne 0 ]
then
    exit 255
fi


# Get the region defined in the current configuration (default to us-west-2 if none defined)
region=$(aws configure get region)
region=${region:-eu-west-2}

ecr_repo="${account}.dkr.ecr.${region}.amazonaws.com"
remote_tag="quay.io/pathogen-genomics-cymru/${image}:${version}"
ecr_tag="${ecr_repo}/tb-pipeline/${image}:${version}"


echo "AWS Region      : ${region}"
echo "Source image tag: $remote_tag"
echo "Target image tag: $ecr_tag"


# If the repository doesn't exist in ECR, create it.
aws ecr describe-repositories --repository-names "tb-pipeline/${image}" > /dev/null 2>&1

if [ $? -ne 0 ]
then
    echo "The repository with name tb-pipeline/${image} does not exist in the registry ${ecr_repo}. Creating repository."
    aws ecr create-repository --repository-name "tb-pipeline/${image}" 
# > /dev/null
fi

# Get the login command from ECR and execute it directly
aws ecr get-login-password --region "${region}" | docker login --username AWS --password-stdin "${account}".dkr.ecr."${region}".amazonaws.com

# Build the docker image locally with the image name and then push it to ECR
# with the full name.

docker pull ${remote_tag}
docker tag ${remote_tag} ${ecr_tag}

docker push ${ecr_tag}
