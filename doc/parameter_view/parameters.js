/*
  Copyright (C) 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


demangleStrings();
reorderList();


// Grab text from the search input field and only show parameters that contain the text:
function filter_text()
{
  var needle = document.getElementById("search").value.toUpperCase();
  var list = document.getElementById("ParameterList");
  var scroll = true;
  collapseAllSubsections();

  var items = document.getElementsByClassName("parameter");
  for (var i=0;i<items.length; ++i) {
    var haystack = items[i].innerText.toUpperCase();
    if (haystack.includes(needle)) {
      items[i].style.display = "block";

      p = items[i];
      while (p.parentNode) {
	p=p.parentNode;
	if (p.className=="content") {
	  var temp = [];
	  temp[0]=p.parentNode.children[0];
	  expand(temp);
	}
      }

      if (scroll) {
	items[i].scrollIntoView(false);
	scroll = false;
      }
    } else {
      items[i].style.display = "none";
    }
  }

  if (needle.length==0) {
    // If the user typed nothing, collapse all
    collapseAllSubsections();
  }

  document.getElementById("link").href="parameters.xml?s=" + encodeURIComponent(document.getElementById("search").value);
}

// Trigger a filter when the user presses return:
function search_input(event) {
  if (event.keyCode == 13)
    filter_text();
}

document.getElementById("search").addEventListener('change', filter_text);
document.getElementById("search").addEventListener('input', search_input);


function demangleStrings() {
  mangledStrings = document.getElementsByClassName("mangled");
  var j;
  for (j = 0; j < mangledStrings.length; j++) {
    mangledStrings[j].innerHTML = mangledStrings[j].innerHTML.replace(/_20/g," ");
    mangledStrings[j].innerHTML = mangledStrings[j].innerHTML.replace(/_2d/g,"-");
  }
  return;
}

function expand(collection) {
  var i;
  for (i = 0; i < collection.length; i++) {
    collection[i].classList.add("active");
    var content = collection[i].nextElementSibling;
    content.style.display = "block";
  }
}

function collapse(collection) {
  var i;
  for (i = 0; i < collection.length; i++) {
    collection[i].classList.remove("active");
    var content = collection[i].nextElementSibling;
    content.style.display = "none";
  }
}

function expandAll() {
  var coll = document.getElementsByClassName("collapsible");
  expand(coll)
}

function collapseAll() {
  var coll = document.getElementsByClassName("collapsible");
  collapse(coll)
}

function expandAllSubsections() {
  var coll = document.getElementsByClassName("subsection");
  expand(coll)
}

function collapseAllSubsections() {
  var coll = document.getElementsByClassName("subsection");
  collapse(coll)
}

function getParameterByName(name, url) {
  if (!url) url = window.location.href;
  name = name.replace(/[\[\]]/g, "\\$&");
  var regex = new RegExp("[?&]" + name + "(=([^&#]*)|&|#|$)"),
      results = regex.exec(url);
  if (!results) return null;
  if (!results[2]) return '';
  return decodeURIComponent(results[2].replace(/\+/g, " "));
}

function autorun() {
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

  // if we get started with parameters.xml?s=bla, perform search
  var search = getParameterByName("s");
  if (search) {
    console.log(search)
    document.getElementById("search").value=search;
    filter_text();
  }
}

function sortTopNodes(ClassType) {
  var list, i, switching, shouldSwitch;
  list = document.getElementById("ParameterList");
  switching = true;

  /* Make a loop that will continue until
     no switching has been done: */
  while (switching) {
    // Start by saying: no switching is done:
    switching = false;
    children = list.children;
    // parameters = children.getElementbyClassName("parameter");

    // Loop through all list items:
    for (i = 0; i < (children.length - 1); i++) {
      if (children[i].children[0].classList.contains(ClassType) && children[i+1].children[0].classList.contains(ClassType)) {
	// Start by saying there should be no switching:
	shouldSwitch = false;
	/* Check if the next item should
	   switch place with the current item: */
	if (children[i].children[0].innerHTML.toLowerCase() > children[i + 1].children[0].innerHTML.toLowerCase()) {
	  /* If next item is alphabetically lower than current item,
	     mark as a switch and break the loop: */
	  shouldSwitch = true;
	  break;
	}
      }
    }

    if (shouldSwitch) {
      /* If a switch has been marked, make the switch
	 and mark the switch as done: */
      children[i].parentNode.insertBefore(children[i + 1], children[i]);
      switching = true;
    }
  }
}

function reorderList() {
  var list, i;
  list = document.getElementById("ParameterList");

  /* Move parameters to the front */
  children = list.children;

  for (i = 0; i < children.length; i++) {
    if (children[i].children[0].classList.contains("parameter")) {
      children[i].parentNode.insertBefore(children[i], children[0]);
    }
  }

  sortTopNodes("parameter");
  sortTopNodes("subsection");
}



window.onload = autorun;
