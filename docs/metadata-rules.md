## BASE
### **`MUST`** /Library: Library has format version
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003186 | False | library format version | False |

### **`MUST`** /Library: Library has name
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003188 | False | library name | False |

### **`MUST`** /Library/Spectrum: Spectrum has unique key
Combination logic: `OR`

| accession | allow_children | name | repeatable | value |
|-----|-----|-----|-----|-----|
| MS:1003237 | False | library spectrum key | False | value_is_unique |

### **`MUST`** /Library/Spectrum: Spectrum has index
Combination logic: `OR`

| accession | allow_children | name | repeatable | value |
|-----|-----|-----|-----|-----|
| MS:1003062 | False | library spectrum index | False | value_is_unique |

### **`MUST`** /Library/Spectrum: Spectrum has index
Combination logic: `OR`

| accession | allow_children | name | repeatable | value |
|-----|-----|-----|-----|-----|
| MS:1003062 | False | library spectrum index | False | value_is_unique |

### **`SHOULD`** /Library/Spectrum: Spectrum has precursor charge
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000041 | False | charge state | False |
| MS:1000633 | False | possible charge state | False |

### **`SHOULD`** /Library/Spectrum: Spectrum has precursor mz
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000744 | False | selected ion m/z | False |
| MS:1003208 | False | experimental precursor monoisotopic m/z | False |

### **`SHOULD`** /Library/Spectrum: Spectrum has aggregation
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003065 | True | spectrum aggregation type | False |

### **`SHOULD`** /Library/Spectrum/Analyte: Analyte has any mass
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1001117 | True | theoretical mass | False |
| MS:1000224 | True | molecular mass | False |

## CONSENSUS
### **`SHOULD`** /Library/Spectrum: Spectrum has replicates used
Combination logic: `AND`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003069 | False | number of replicate spectra available | False |
| MS:1003070 | False | number of replicate spectra used | False |

## GOLD
### **`SHOULD`** /Library: Library has contact
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000586 | False | contact name | False |
| MS:1000587 | False | contact address | False |
| MS:1000588 | False | contact URL | False |
| MS:1000589 | False | contact email | False |
| MS:1000590 | False | contact affilation | False |

### **`SHOULD`** /Library: Library has reference
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1002866 | False | Reference | False |

## PEPTIDE
### **`MUST`** /Library/Spectrum/Analyte: Analyte has peptide seq
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000889 | False | proforma peptidoform sequence | False |
| MS:1000888 | False | stripped peptide sequence | False |
| MS:1003270 | False | proforma peptidoform ion notation | False |

## SILVER
### **`SHOULD`** /Library: Library has identifier
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003187 | False | library identifier | False |

### **`SHOULD`** /Library/Spectrum: Spectrum has origin type
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003072 | True | spectrum origin type | False |

### **`SHOULD`** /Library/Spectrum: Spectrum has dissociation
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000044 | True | dissociation method | False |

## SINGLE
### **`SHOULD`** /Library/Spectrum: Spectrum has source file
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1000577 | False | source file | False |

### **`SHOULD`** /Library/Spectrum: Spectrum has scan identifier
Combination logic: `OR`

| accession | allow_children | name | repeatable |
|-----|-----|-----|-----|
| MS:1003057 | False | scan number | False |

