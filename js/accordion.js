$(document).ready(function () {
  $("body").ready(function () {
    $(".accordion > .foldable > .fold-content").css("display", "none");
    $(".accordion > .foldable.fold-default > .fold-label").addClass("open");
    $(".accordion > .foldable.fold-default > .fold-content").slideDown(350);
  })
  $(".accordion > .foldable > .fold-label").on("click", function (e) {
    if($(this).parent().has(".accordion")) {
      e.preventDefault();
    }
    if(!$(this).hasClass("open")) {
      $(".accordion > .foldable > .fold-content").slideUp(350);
      $(".accordion > .foldable > .fold-label").removeClass("open");

      $(this).addClass("open");
      $(this).next(".fold-content").slideDown(350);
    }
    
    else if($(this).hasClass("open")) {
      $(this).next(".fold-content").slideUp(350);
      $(this).removeClass("open");
    }
  });
});
