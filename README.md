# fgsea operator

##### Description
`fgsea` operator performs a fast gene set enrichment analysis.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, input data usually a ranking statistics, per cell

Output relations|.
---|---
`NES`         | numeric, fgsea of the input data
`padj`        | numeric, adjusted p-value

##### Details

The operator takes all the values of a cell and returns the value which is the fgsea. The computation is done per cell. There is one value returned for each of the input cell.

#### References


##### See Also


#### Examples
