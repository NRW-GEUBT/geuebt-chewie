# Guide for developers

## Changing databases

The easiest way to use alternative databases is to download them manually and
provide them as paths in the configuarion.

If you wish to deploy the module with other databases, you can also fork the
repository and modify the post-deploy scripts for each environement
(`<ENV NAME>.post.deploy.sh`).
Simply modifying the URL vairable should be enough in most cases.

## Data validation

Data validation is performed through JSON schema validations.
The schemas used for validation are in the `workflow/schema` directory.
Validation is perfromed with the [jsonschema](https://python-jsonschema.readthedocs.io/en/latest/#)
python package and follows the [2020-12](https://json-schema.org/draft/2020-12/release-notes.html) specification.