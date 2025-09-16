## Library for `IN.DAT` Files

> `process/io/in_dat.py`

A set of Python classes to read, modify and write an `IN.DAT` file. 

### `class InDat`

Class `InDat` for handling `IN.DAT` data. It handles

- Reading IN.DAT files
- Writing IN.DAT files
- Storing information in dictionary for use in other codes
- Alterations to IN.DAT

To open a PROCESS input file and store the information in a Python class use:

```python
    i = InDat(filename="IN.DAT")
```

To get values of input file entries from the object

```python
    i.data["ixc"].value
    i.number_of_constraints
    i.number_of_itvars
    i.data["fimp"].value
    i.data["i_pf_location"].value
```

To add/remove constraints, iteration variables

```Python
    i.remove_constraint_equation(2)
    i.add_constraint_equation("3")
    i.add_constraint_equation(2)
    i.add_iteration_variable(103)
    i.add_iteration_variable("2")
    i.remove_iteration_variable(2)
    i.remove_iteration_variable("3")
    # Add bound will change the bound value if it already exists
    i.add_bound(103, "upper", 5.0)
    i.remove_bound(2, "upper")
```

To change the value or add/remove a regular input

```python
# Add parameter will change the parameter value if it already exists
i.add_parameter("blnktthdsd", 0.5)
i.add_parameter("iavail", 1)
i.remove_parameter("blnkithsddd")
i.remove_parameter("dr_blkt_inboard")
i.add_parameter("sweep", [3.0, 3.0])
```

To write the modified `IN.DAT` to file

```python
i.write_in_dat()
i.write_in_dat(output_filename="new_IN.DAT")
```

## Library for `MFILE.DAT` Files

> `process/io/mfile.py`


A set of Python classes to read and extract data from the `MFILE.DAT`.

### `class MFile`

To open an MFILE.DAT as a Python object:

```python
m = MFile(filename="MFILE.DAT")
```

As a PROCESS run may contain a scan, access to the MFILE data is determined 
by scan number. To find the number of scans use

```python
m.data["rmajor"].get_number_of_scans()
```

To get an individual scan for a variable use

```python
m.data["rmajor"].get_scan(2)
```

If the file contains no scan, but a single data point then one can use the above 
with the argument of `-1` like below

```python
m.data["rmajor"].get_scan(-1)
```

To get all of the scans as a list (if there is a single point it returns list 
length of 1) use

```python
m.data["rmajor"].get_scans()
```
