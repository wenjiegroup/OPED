function exam_seq_fun (obj) {
    switch(obj) {
        case 1:
            document.getElementById("TEXT_SEQ").value = "CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTGGCGCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC";
            break;
        case 2:
            document.getElementById("TEXT_SEQ").value = "CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGCGCTG(+ATT)GCGCGAGGCCGCCTGGCAACTCTGCGACTACTACCTGCC";
            break;
        case 3:
            document.getElementById("TEXT_SEQ").value = "CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCGCGGCCGCGTCAAGCCGGGGGGGC(-GCTGGCGCGA)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC";
            break;
    }
}

function exam_pos_fun (obj) {
    document.getElementById("Chromosome").value = "chr1";
    document.getElementById("Position").value = "943995";
    switch(obj) {
        case 1:
            document.getElementById("Pattern").value = "/T";
            break;
        case 2:
            document.getElementById("Pattern").value = "+GTATT";
            break;
        case 3:
            document.getElementById("Pattern").value = "-GAGAACTCGG";
            break;
    }
}


function exam_var_fun (obj) {
    switch(obj) {
        case 1:
            document.getElementById("queryType").value = "AlleleID";
            document.getElementById("queryItem").value = "929884";
            break;
        case 2:
            document.getElementById("queryType").value = "GeneID";
            document.getElementById("queryItem").value = "148398";
            break;
        case 3:
            document.getElementById("queryType").value = "GeneSymbol";
            document.getElementById("queryItem").value = "SAMD11";
            break;
        case 4:
            document.getElementById("queryType").value = "HGNC_ID";
            document.getElementById("queryItem").value = "28706";
            break;
    }
}


