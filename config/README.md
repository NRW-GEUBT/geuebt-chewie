# Workflow configuration

## Configuration

To configure the workflow, modifiy `config/config.yaml` according
to your needs.

For specific details about the parameters
check the documentation.

## Additional files

Additional files must be provided as a tab-delimited text data, with strict requirements:

### Sample sheet

| Field name             | Value                               |
| ---------------------- | ----------------------------------- |
| sample                 | Sample name                         |
| assembly               | Path to the assembly file           |


### External clusters

| Field name             | Value                               |
| ---------------------- | ----------------------------------- |
| sample                 | Sample name                         |
| cluster_name           | Name of the cluster in the externa database.  |

The cluster naming should follow the structure `<PREFIX>-<Num>.<Num>`. For example
LIS-125.32, meaning Listeria cluster number 125, subcluster 32.
In the main cluster file, only provide the names trucncated to the cluster number (LIS-125),
and the full name in the subcluster file (LIS-125.32).

### External profiles

| Field name             | Value                               |
| ---------------------- | ----------------------------------- |
| #FILE                  | Sample name                         |
| <allele file name>     | CRC32 hash of the corresponding allele  |

Truncated example:

```plaintext
#FILE    lmo0001.fasta    lmo0002.fasta    lmo0003.fasta    lmo0005.fasta
DB-LI00339-0    3453202319    138852938    2562060452    1874597132    2322454825
DB-LI00345-0    3453202319    138852938    874063884    1874597132    2322454825
```

### External timestamps

| Field name             | Value                               |
| ---------------------- | ----------------------------------- |
| sample                 | Sample name                         |
| date                   | Sample timestamp (ISO8601, e.g. 2023-04-30) |

### External statistics

| Field name             | Value                               |
| ---------------------- | ----------------------------------- |
| sample                 | Sample name                         |
| EXC                   | Same named parameter from CHEWBBACCA results |
| INT                   | Same named parameter from CHEWBBACCA results |
| LNF                   | Same named parameter from CHEWBBACCA results |
| PLOT                   | Same named parameter from CHEWBBACCA results |
| NIPH                   | Same named parameter from CHEWBBACCA results |
| ALM                   | Same named parameter from CHEWBBACCA results |
| ASM                   | Same named parameter from CHEWBBACCA results |
