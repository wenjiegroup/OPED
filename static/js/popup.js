// get pop-up window
// https://c.runoob.com/codedemo/2041/
var modal = document.getElementById('myModal');
 
// open window button
// var btn = document.getElementById("myBtn");
 
// get <span>，to close window
var span = document.querySelector('.close');
 
// open window
// btn.onclick = function() {
//     modal.style.display = "block";
// }

// click other place, close window
window.onclick = function(event) {
    if (event.target == modal) {
        modal.style.display = "none";
    }
}

// click <span> (x), close window
span.onclick = function(obj) {
    modal.style.display = "none";
}