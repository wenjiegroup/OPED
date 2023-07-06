var igvDiv=document.getElementById("igv_div")
var options = {
            reference: {
                fastaURL: "test.fasta"
            }
        };

        igv.createBrowser(igvDiv, options)
            .then(function(browser) {
                console.log("Created IGV browser");
            })
