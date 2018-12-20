demangleStrings();

function demangleStrings() {
mangledStrings = document.getElementsByClassName("mangled");
var j;
for (j = 0; j < mangledStrings.length; j++) {
mangledStrings[j].innerHTML = mangledStrings[j].innerHTML.replace(/_20/g," ");
mangledStrings[j].innerHTML = mangledStrings[j].innerHTML.replace(/_2d/g,"-");
}
return;
}


var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
}
