    function toggleMe(btn,a){
      var e=document.getElementById(a);
      if(!e)return true;
      if(e.style.display=="none"){
        e.style.display="block";
        btn.value = "-";
      }
      else{
        e.style.display="none";
        btn.value = "+";
      }
       return true;
    }
    function showbtn(){
      var btns = document.querySelectorAll(".btn");
      if(!btns) alert('????');
      for(var btn of btns){
        btn.style.display="inline";
        btn.style.fontFamily="monospace";
        btn.style.color="Maroon";
        btn.style.margin="0px 0px 0px 0px";
        btn.style.paddingRight="3px";
        btn.style.paddingLeft="3px";
/*  border:0px none Maroon


  padding-left: 5px;
  text-decoration: none;*/

      }

    }
