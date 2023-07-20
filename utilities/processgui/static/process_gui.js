function highlight(obj) {
    var refname = "ref_" + obj.name
    ref = $(obj).parent().parent().find("[name='" + refname + "']")
    if(obj.value != ref.val()){
        obj.style.color = "#ff0000";
        ref.css("color", "#ff0000");
        }
    else {
        obj.style.color = "#000000";
        ref.css("color", "#000000");
        }
    return true;
    }


//set the value of all objects with the same name
//as the given object to the value of the given object
function syncval(obj) {
    var name=obj.name
    var value=obj.value
    var sel = $("[name='"+name+"']")
    sel.val(value);
    $.each(sel, function() {highlight(this)})
    }

function showall() {
    $( ".accordion" ).accordion({
    active: 0,
    });
}


function hideall() {
    $( ".accordion" ).accordion({
        active: false
    });
}
