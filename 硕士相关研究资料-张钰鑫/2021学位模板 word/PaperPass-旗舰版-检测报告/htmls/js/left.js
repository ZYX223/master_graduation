(function(System, $) {
    $(function() {



        parent.parent.postMessage("word", "*");
        var storage = window.localStorage;
        var LineState = System.operateState.getLineState(Report.report_id);
        var TypeOf = typeof isHasDetailReport;
        var ProTypeOf = typeof isHasProofReport;
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            var DetailNum = System.operateState.getDetailNum(Report.report_id);
            var Detailoffset = System.operateState.getDetailoffset(Report.report_id);
            var DetailSize = System.operateState.getDetailSize(Report.report_id);
            var DetailTopNum = System.operateState.getDetailTopNum(Report.report_id);
            var DetailPercent = System.operateState.getDetailPercent(Report.report_id);
            var ZoomBtnState = System.operateState.getZoomBtnState(Report.report_id);
            var FoldState = System.operateState.getFoldState(Report.report_id);
        } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
            var ProofNum = System.operateState.getProofNum(Report.report_id);
            var ProofOffset = System.operateState.getProofOffset(Report.report_id);
            var ProofSize = System.operateState.getProofSize(Report.report_id);
            var ProofTopNum = System.operateState.getProofTopNum(Report.report_id);
            var ProofPercent = System.operateState.getProofPercent(Report.report_id);
            var ProZoomBtnState = System.operateState.getProZoomBtnState(Report.report_id);
            var FoldState2 = System.operateState.getFoldState2(Report.report_id);
        } else {
            var OriNum = System.operateState.getOriNum(Report.report_id);
            var Orioffset = System.operateState.getOrioffset(Report.report_id);
            var OriSize = System.operateState.getOriSize(Report.report_id);
            var OriTopNum = System.operateState.getOriTopNum(Report.report_id);
            var OrilPercent = System.operateState.getOrilPercent(Report.report_id);
            var OriZoomBtnState = System.operateState.getOriZoomBtnState(Report.report_id);
        };

        var str = "";
        str += `<span class="pageNum"></span><span class="partition"></span><span class="totalPage"></span>`;
        $(".top_Num").html(str);

        function FunPage(PageNum, offsetH) {
            if (PageNum == "null" || PageNum == "") {
                $(".pageNum").text(1);
                page1 = $(".input_pageNum").val(1);
            } else {
                $(".pageNum").text(PageNum);
                $(".input_pageNum").val(PageNum);
                $("#toall_Num").attr("href", "#" + PageNum);
                $("#body").scrollTop(offsetH);
            };
        };
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            parentWidth(DetailSize ? DetailSize : 1);
            AddH(DetailSize ? DetailSize : 1);
            FunPage(DetailNum, Detailoffset);
            size = (DetailSize ? DetailSize : 1);

            if (ZoomBtnState == "" || ZoomBtnState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#323639");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
                $("#menuButton").css("background", "#1b1e1f");
            };




        } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
            parentWidth(ProofSize ? ProofSize : 1);
            AddH(ProofSize ? ProofSize : 1);
            FunPage(ProofNum, ProofOffset);
            size = (ProofSize ? ProofSize : 1);
            if (ProZoomBtnState == "" || ProZoomBtnState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#323639");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
                $("#menuButton").css("background", "#1b1e1f");
            };
        } else {
            parentWidth(OriSize ? OriSize : 1);
            AddH(OriSize ? OriSize : 1);
            FunPage(OriNum, Orioffset);
            size = (OriSize ? OriSize : 1);

            if (OriZoomBtnState == "" || OriZoomBtnState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#323639");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
                $("#menuButton").css("background", "#1b1e1f");
            };




        };
        $(".partition").text("/");
        $(".top_Num").css({
            "background": "#000000",
            "height": "30px",
            "line-height": "30px",
            "padding": "0 10px",
        });
        $(".top_Num").show();

        $(".Switch").click(function() {
            $(".layui-layer-shade").css("display", "block");
            $(".layui-layer").css("display", "block");
            $(".layui-layer-title").html("切换模式");
            $(".layui-layer-content").html("确认切换至“纯文字模式”吗？");
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                $(".detail_btn").data('type', 'switch_text');
            }
            if (ProTypeOf == "boolean" && isHasProofReport == true) {
                $(".detail_btn").data('type', 'proof_switch_text');
            }
        });

        $(".detail_btn").click(function() {
            type = $(this).data('type');
            if (type == 'switch_text') {
                $(".layui-layer-btn0").attr("href", $(this).data('txtLink'));
                parent.parent.postMessage("report_text", "*");
                parent.postMessage("word_right", "*");
                System.operateState.saveInitPage(Report.report_id, type);
            }
            if (type == 'proof_switch_text') {
                $(".layui-layer-btn0").attr("href", $(this).data('txtLink'));
                parent.postMessage("proof_right", "*");
                parent.parent.postMessage("proof_text", "*");
                System.operateState.saveProofInitPage(Report.report_id, type);
            }
        });
        $(".layui-layer-setwin").click(function() {
            $(".layui-layer-shade").css("display", "none");
            $(".layui-layer").css("display", "none");
        });
        $(".layui-layer-btn1,.layui-layer-btn0 ").click(function() {
            $(".layui-layer-shade").css("display", "none");
            $(".layui-layer").css("display", "none");
        });

		$(".subline").click(function () {
			$(".layui-layer-shade").css("display", "block");
			$(".layui-layer").css("display", "block");
			$(".layui-layer-title").html("辅助线");
			curHref = $(".detail_btn").attr("href");
			if (curHref != '#') {
				$(".detail_btn").data('txtLink', $(".detail_btn").attr("href"));
				$(".layui-layer-btn0").attr("href", "#");
			};
			$(".detail_btn").data('type', 'subline');
			if ($(".red").parent().eq(0).hasClass('red_a') || $(".orange").parent().eq(0).hasClass('orange_a')) {
				$(".layui-layer-content").html("确认“取消辅助线”吗？");
			} else {
				$(".layui-layer-content").html("确认“使用辅助线”吗？");
			};
		});
		function FuncLine(Line) {
			if (Line == "true" || Line == "") {
				$(".red").parent().removeClass('red_a');
				$(".orange").parent().removeClass('orange_a');
				$(".subline").css("background-position", "-106px -162px");
			} else {
				$(".red").parent().addClass('red_a');
				$(".orange").parent().addClass('orange_a');
				$(".subline").css("background-position", "-146px -159px");
			};
		};
		$('body').on('click', '.layui-layer-btn0', function () {
			if ($(".layui-layer-title").html() !== '辅助线') {
				return;
			};
			var has_red = $("span").hasClass("red");
			var has_orange = $("span").hasClass("orange");
			if (has_red == false && has_orange == false) {
				FuncLine("");
			} else {
				var has_Class = $(".red").parent().eq(0).hasClass('red_a') || $(".orange").parent().eq(0).hasClass('orange_a');
				var str_has = String(has_Class);
				System.operateState.saveLineState(Report.report_id, str_has);
				FuncLine(str_has);
			};
		});
		FuncLine(LineState);

        var totalPage = $(".upkmzjakmc").length;
        $('.totalPage').text(totalPage);

        var nowPage = 1;

        function scrollFunc(size) {
            var offsetH = $("body").scrollTop();
            var scrollH = -15;
            $('.upkmzjakmc').each(function(k, v) {
                var PageH = $(v)[0].getBoundingClientRect().height + 15;
                if (scrollH < offsetH && offsetH <= (scrollH + PageH)) {
                    nowPage = k + 1;
                };
                scrollH += PageH;
            });
            $('.input_pageNum').val(nowPage);
            $('.pageNum').text(nowPage);
            $("#toall_Num").attr("href", "#" + nowPage);
        };

        function parentWidth(size) {
            var _Size = Math.round(size * 100) / 100;
            percent = Math.round(_Size * 100);
            if (percent == "") {
                $(".percent_num").text(100);
            } else {
                $(".percent_num").text(percent);
            };

            $(".upkmzjakmc").css("transform", "scale(" + _Size + ")");
            if (_Size < 1) {
                $(".upkmzjakmc").css("transform-origin", "left top");
            } else {
                $(".upkmzjakmc").css("transform-origin", "center top");
            }
            if (nowPage != undefined) {
                $(".jump_page").attr("href", "#" + nowPage);
            };
            $('.upkmzjakmc').each(function(k, v) {
                var pageW = ($(this)[0].getBoundingClientRect().width) + 30;
                var pageH = ($(this)[0].getBoundingClientRect().height) + 15;
                $(this).parent(".jcjymywwzd").css("width", pageW);
                $(this).parent(".jcjymywwzd").css("height", pageH);
            });
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                $("body").off("scroll").on("scroll", function(e) {
                    scrollFunc(_Size);
                    var offsetH = $("body").scrollTop();
                    var Page = nowPage;
                    var obj = {
                        PageNum1: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.saveDetailNum(Report.report_id, Page);
                    System.operateState.saveDetailoffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.saveDetailSize(Report.report_id, PageSize);
                System.operateState.saveDetailPercent(Report.report_id, percent);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                $("body").off("scroll").on("scroll", function(e) {
                    scrollFunc(_Size);
                    var offsetH = $("body").scrollTop();


                    var Page = nowPage;
                    var obj = {
                        PageNum5: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.saveProofNum(Report.report_id, Page);
                    System.operateState.saveProofOffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.saveProofSize(Report.report_id, PageSize);
                System.operateState.saveProofPercent(Report.report_id, percent);
            } else {
                $("body").off("scroll").on("scroll", function(e) {
                    scrollFunc(_Size);
                    var offsetH = $("body").scrollTop();
                    var Page = nowPage;
                    var obj = {
                        PageNum2: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.saveOriNum(Report.report_id, Page);
                    System.operateState.saveOrioffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.saveOriSize(Report.report_id, PageSize);
                System.operateState.saveOrilPercent(Report.report_id, percent);
            };
        };

        function TopNum() {
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                var TopNum = $(".top_Num").css("right");
                System.operateState.saveDetailTopNum(Report.report_id, TopNum);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                var TopNum = $(".top_Num").css("right");
                System.operateState.saveProofTopNum(Report.report_id, TopNum);
            } else {
                var TopNum = $(".top_Num").css("right");
                System.operateState.saveOriTopNum(Report.report_id, TopNum);
            };
        };

        function HideShade() {
            setTimeout(function() {
                $(".zoom-shade").hide();
            }, 200);
        };
        var percent;

        $("#originalButton").on("click", function() {
            $(".zoom-shade").show();
            size = 1;
            _size = size;
            parentWidth(_size);
            AddH(_size);
            TopNum();
            HideShade();
        });

        $("#fitPageButton").click(function() {
            $(".zoom-shade").show();
            var fitPage_h = $("body")[0].getBoundingClientRect().height;
            size = (fitPage_h - 20) / minHeight;
            _size = size;
            parentWidth(_size);
            AddH(_size);
            TopNum();
            HideShade();
        });

        $("#fitWidthButton").click(function() {
            $(".zoom-shade").show();
            var fitWidth_w = $("body")[0].getBoundingClientRect().width;
            size = (fitWidth_w - 163) / minWidth;
            if (size > 2.0) {
                size = 2.0;
            };
            _size = size;
            parentWidth(_size);
            AddH(_size);
            TopNum();
            HideShade();
        });

        function roundFun(value, n) {
            return Math.round(value * Math.pow(10, n)) / Math.pow(10, n);
        };

        $("#zoomOutButton").click(function() {
            $(".zoom-shade").show();
            if (percent % 10 == 0 && size <= 2.0) {
                size += 0.1;
            };
            if (percent % 10 >= 5 && size <= 2.0) {
                size = roundFun(percent / 100, 1);
            };
            if (percent % 10 < 5 && size <= 2.0) {
                size = roundFun(percent / 100, 1) + 0.1;
            };
            if (size > 2.0) {
                size = 2.0;
            };
            _size = size;
            parentWidth(_size);
            AddH(_size);
            TopNum();
            HideShade();
        });

        $("#zoomInButton").click(function() {
            $(".zoom-shade").show();
            if (percent % 10 == 0 && size >= 0.1) {
                size -= 0.1;
            } else if (percent % 10 >= 5 && size >= 0.1) {
                size = roundFun(percent / 100, 1) - 0.1;
            } else if (percent % 10 < 5 && size >= 0.1) {
                size = roundFun(percent / 100, 1);
            };
            if (size < 0.1) {
                size = 0.1;
            };
            _size = size;
            parentWidth(_size);
            AddH(_size);
            TopNum();
            HideShade();
        });

        $(".allButton").mousemove(function() {
            $(this).css("background", "#1b1e1f");
            $(this).children(".bot_tip").css("display", "block");
        }).mouseout(function() {
            $(this).css("background", "#323639");
            $(this).children(".bot_tip").css("display", "none");
        });

        $("#CollagenButton").click(function(e) {
            var listState = $(".percent_list").css("display");
            if (listState == "block") {
                $(".percent_list").hide();
                $(this).children(".percent_tip").css("display", "block");
            } else {
                $(".percent_list").show();
                $("#CollagenButton").css("background", "#1b1e1f");
            };
            $(this).children(".percent_tip").css("display", "none");
            e.stopPropagation();
        }).mousemove(function() {
            $("#CollagenButton").css("background", "#1b1e1f");
            if ($(".percent_list").css("display") == "block") {
                $(this).children(".percent_tip").css("display", "none");
            } else {
                $(this).children(".percent_tip").css("display", "block");
            };
        }).mouseout(function() {
            if ($(".percent_list").css("display") == "block") {
                $("#CollagenButton").css("background", "#1b1e1f");
            } else {
                $("#CollagenButton").css("background", "#323639");
            };
            $(this).children(".percent_tip").css("display", "none");
        });
        $(document).click(function() {
            $(".percent_list").hide();
            $("#CollagenButton").css("background", "#323639");
        });
        $(".zoom-shade").click(function(e) {
            var listState = $(".percent_list").css("display");
            if (listState == "block") {
                $(".percent_list").show();
            } else {
                $(".percent_list").hide();
            };
            e.stopPropagation();
        });

        $(".percent_li").click(function(e) {
            $(".zoom-shade").show();
            var liTxt = $(this).children(".li_text").text();
            size = liTxt / 100;
            _size = size;
            parentWidth(_size);
            AddH(_size);
            HideShade();
            TopNum();
            e.stopPropagation();
        }).mousemove(function() {
            $(this).css("background", "#1b1e1f")
        }).mouseout(function() {
            $(this).css("background", "#323639");
        });;

        $("#menuButton").click(function() {


            var showState = $(".allButton").css("display");
            $(".percent_list").hide();
            $("#CollagenButton").css("background", "#323639");
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                var showState = $(".allButton").css("display");
                System.operateState.saveZoomBtnState(Report.report_id, showState);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                var showState = $(".allButton").css("display");
                System.operateState.saveProZoomBtnState(Report.report_id, showState);
            } else {
                var showState = $(".allButton").css("display");
                System.operateState.saveOriZoomBtnState(Report.report_id, showState);
            };
            if (showState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#323639");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
                $("#menuButton").css("background", "#1b1e1f");
            };








        }).mousemove(function() {
            $("#menuButton").css("background", "#1b1e1f")
            $(this).children(".menu_tip").css("display", "block");
        }).mouseout(function() {
            $(this).children(".menu_tip").css("display", "none");
            if ($(".allButton").css("display") == "block") {
                $("#menuButton").css("background", "#1b1e1f");
            } else {
                $("#menuButton").css("background", "#323639");
            };
        });
        $(".Tip-page .layui-layer-setwin").click(function() {
            $(".Tip-shade").hide();
            $(".Tip-page").hide();
        });
        $(".Tip-page .layui-layer-btn").click(function() {
            $(".Tip-shade").hide();
            $(".Tip-page").hide();
        });


        $(".bottom").click(function() {
            $(".bottom_a").attr("href", "#" + (nowPage + 1));
            $("#toall_Num").attr("href", "#" + (nowPage + 1));
        });
        $(".top").click(function() {
            $(".top_a").attr("href", "#" + (nowPage - 1));
            $("#toall_Num").attr("href", "#" + (nowPage - 1));
        });

        $(".input_pageNum").focus(function() {
            $(".input_pageNum").css("background-color", "#424649");
            $(".input_pageNum").select();
        }).blur(function() {
            $(".input_pageNum").css("background-color", "transparent");
            page1 = $(".input_pageNum").val();
            if (page1 <= totalPage && /^[1-9]\d*$/.test(page1)) {
                $("#toall_Num").attr("href", "#" + page1);
            } else {
                $(".input_pageNum").val(nowPage);
            };
        }).mousemove(function() {
            $(".input_pageNum").css("background-color", "#424649");
        }).mouseout(function() {
            if ($('.input_pageNum').is(':focus')) {
                $(".input_pageNum").css("background-color", "#424649");
            } else {
                $(".input_pageNum").css("background-color", "transparent");
            };
        }).keyup(function() {
            if (event.keyCode == 13) {
                page1 = $(".input_pageNum").val();
                if (page1 <= totalPage && /^[1-9]\d*$/.test(page1)) {
                    $("#toall_Num").attr("href", "#" + page1);
                    window.location.hash = "#" + page1;
                } else {
                    $(".input_pageNum").val(nowPage);
                    $(".input_pageNum").blur();
                };
            };
        });

        $(".toall_div").mousemove(function() {
            $(this).css({
                "background-color": "#424649",
                "border-radius": "2px"
            });
            $(this).children(".tip").css("display", "block");
        }).mouseout(function() {
            $(this).css("background-color", "transparent");
            $(this).children(".tip").css("display", "none");
        });

        var w1 = screen.width;
        var h1 = screen.height;
        var w = document.documentElement.clientWidth;
        var h = document.documentElement.clientHeight;
        var left = (w - minWidth * size) / 2;
        var leftWidth = left - 75;
        var rightWidth = left + 15;
        var shadeLeft = (w - 500) / 2;
        var shadePercent = (w - 500) / 2;
        var shadeTop = (h1 - 273) / 2;
        $(".Shadow ").css("width", w - 17);
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            if (minWidth * size >= (w - 105)) {
                $(".top_Num").css("right", 77);
            } else {
                if (DetailTopNum == "") {
                    $(".top_Num").css("right", rightWidth);
                } else {
                    $(".top_Num").css("right", DetailTopNum);
                };
            };
        } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
            if (minWidth * size >= (w - 105)) {
                $(".top_Num").css("right", 77);
            } else {
                if (ProofTopNum == "") {
                    $(".top_Num").css("right", rightWidth);
                } else {
                    $(".top_Num").css("right", ProofTopNum);
                };
            };
        } else {
            if (minWidth * size >= (w - 105)) {
                $(".top_Num").css("right", 77);
            } else {
                if (OriTopNum == "") {
                    $(".top_Num").css("right", rightWidth);
                } else {
                    $(".top_Num").css("right", OriTopNum);
                };
            };
        };
        $("#body").height(h);
        $(".overflow_ul").css("max-height", h - 100);
        $(".layui-nav").css("max-width", w * 0.29);
        $(".layui-layer").css({
            "left": shadeLeft,
            "top": "60px"
        });
        $(".Tip-page").css({
            "left": shadePercent,
            "top": "60px"
        });
        var arr = JSON.parse(detail_title_info);
        var str = `<div class="similarTip" style="position:absolute;border:1px solid #000000;background:#ffffff;z-index:9999;font-size:13px;padding:6px 4px;white-space:nowrap;"></div>`;
        var htm = `<div class="similarTip" style="position:absolute;border:1px solid #000000;background:#ffffff;z-index:9999;font-size:13px;line-height:22px;padding:3px 5px;">
	 					<div class="num"></div>
						<div>查看句子详细结果，请展开右侧！</div>
                   </div>`;
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            $("#mbxsgxktna a[target='right']").mouseover(function() {
                var FoldState = System.operateState.getFoldState(Report.report_id);
                if (FoldState == 'none') {
                    $('#mbxsgxktna').append(htm);
                } else {
                    $('#mbxsgxktna').append(str);
                }
            }).mouseout(function() {
                $("div[class='similarTip']").remove();
            }).click(function() {
                var FoldState = System.operateState.getFoldState(Report.report_id);
                if (FoldState == 'none') {
                    return false;
                }
            });
            $("#mbxsgxktna a[target='right']").mousemove(function(eve) {
                var w = document.documentElement.clientWidth;
                var h = document.documentElement.clientHeight;
                var e = eve || window.event;
                var a_Href = $(this).attr("href").replace(/[\|\/]/g, "").replace(/..sentence_detail/, "").replace(/.html/, "");
                $(this).parent().css("z-index", "");
                var FoldState = System.operateState.getFoldState(Report.report_id);
                if (FoldState == 'none') {
                    $("div[class='similarTip'] .num").html('相似度: ' + arr[a_Href].score + '%');
                } else {
                    $("div[class='similarTip']").html('相似度: ' + arr[a_Href].score + '%');
                }
                $("div[class='similarTip']").css({ "left": (e.clientX + 20) + "px", "top": (e.clientY + 20) + "px" });
                if (h - e.clientY < 70) {
                    $("div[class='similarTip']").css({ "left": (e.clientX + 20) + "px", "top": (e.clientY - 50) + "px" });
                };
                if (w - e.clientX < 120) {
                    $("div[class='similarTip']").css({ "left": (e.clientX - 100) + "px" });
                };
            });
        }
        

        function getScrollBarWidth() {
            var outer = document.createElement("div");
            outer.style.visibility = "hidden";
            outer.style.width = "100px";
            outer.style.position = "absolute";
            outer.style.top = "-9999px";
            document.body.appendChild(outer);
            var widthNoScroll = outer.offsetWidth;
            outer.style.overflow = "scroll";
            var inner = document.createElement("div");
            inner.style.width = "100%";
            outer.appendChild(inner);
            var widthWithScroll = inner.offsetWidth;
            outer.parentNode.removeChild(outer);
            scrollBarWidth = widthNoScroll - widthWithScroll;
            return scrollBarWidth;
        }
        var width = getScrollBarWidth();
        $('.btn_fold').css("right", width + 'px');

        $(".btn_fold").click(function() {
            var fold = $('.btn_fold_left').css('display');
            foldFun(fold);
        })
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            foldFun(FoldState);
        }else if(ProTypeOf == "boolean" && isHasProofReport == true){
            foldFun(FoldState2);
        }

        function foldFun(FoldState) {
            if (FoldState == "none") {
                $('.btn_fold_left').css("display", "block");
                $('.btn_fold_right').css("display", "none");
            } else {
                $('.btn_fold_left').css("display", "none");
                $('.btn_fold_right').css("display", "block");
            }
            if (TypeOf == "boolean" && isHasDetailReport == true){
                var obj = { FoldState: FoldState };
                System.operateState.saveFoldState(Report.report_id, FoldState);
            }else if(ProTypeOf == "boolean" && isHasProofReport == true){
                var obj = { FoldState2: FoldState };
                System.operateState.saveFoldState2(Report.report_id, FoldState);
            }
            parent.postMessage(obj, "*");
            parent.parent.postMessage(obj, "*");
        };
        var page_nav_click = "false";
        window.addEventListener('message', function(event) {
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                if (event.data.FoldState) {
                    var _size = System.operateState.getDetailSize(Report.report_id);
                    parentWidth(_size);
                    AddH(_size);
                    TopNum();
                    HideShade();
                }
            }
            if (ProTypeOf == "boolean" && isHasProofReport == true){
            	if(event.data.FoldState2){
            		var _size = System.operateState.getProofSize(Report.report_id);
            		parentWidth(_size);
            		AddH(_size);
            		TopNum();
            		HideShade();
            	}
            }
            if (event.data.EnterId) {
                window.location.hash = "#" + event.data.EnterId;
            }
            if (event.data.page_nav_click) {
                page_nav_click = event.data.page_nav_click;
            }
        });
    });
})(Report, jQuery);
var arrW = new Array();
var arrH = new Array();
$('.upkmzjakmc').each(function(k, v) {
    arrW.push($(this).width());
});
$('.upkmzjakmc').each(function(k, v) {
    arrH.push($(this).height());
});

function getMaxCountWidth(arr) {
    var mp = new Array();
    for (var idx in arr) {
        var item = arr[idx];
        var count = 0;
        if (!mp.hasOwnProperty(item)) {
            mp[item] = 1;
        } else {
            count = mp[item];
            count = count + 1;
            mp[item] = count;
        };
    };
    var maxCount = 0;
    var width = 0;
    for (var key in mp) {
        var val = mp[key];
        if (val >= maxCount) {
            maxCount = val;
            if (width == 0) {
                width = key;
            } else if (key > width) {
                width = key;
            };
        };
    };
    return width;
};

function getMaxCountHeight(arr) {
    var mp = new Array();
    for (var idx in arr) {
        var item = arr[idx];
        var count = 0;
        if (!mp.hasOwnProperty(item)) {
            mp[item] = 1;
        } else {
            count = mp[item];
            count = count + 1;
            mp[item] = count;
        };
    };
    var maxCount = 0;
    var height = 0;
    for (var key in mp) {
        var val = mp[key];
        if (val >= maxCount) {
            maxCount = val;
            if (height == 0) {
                height = key;
            } else if (key < height) {
                height = key;
            };
        };
    };
    return height;
};
var minWidth = getMaxCountWidth(arrW);
var minHeight = getMaxCountHeight(arrH);
var size = 1;
var _size;

function AddH(size) {
    var h = document.documentElement.clientHeight;
    var last_h = $(".upkmzjakmc:last")[0].getBoundingClientRect().height + 30;
    var AddH = h - last_h;
    if (last_h >= h) {
        $("#Add_height").css("height", 0);
    } else {
        $("#Add_height").css("height", AddH);
    };
    var w = document.documentElement.clientWidth;
    var mWidth = minWidth * size;
    var topNumWidth = $(".top_Num").width() + 25;
    var right = (w - mWidth) / 2;
    if (mWidth <= 500) {
        $(".top_Num").css("right", right - topNumWidth);
    } else if (mWidth >= (w - 105)) {
        $(".top_Num").css("right", 77);
    } else {
        $(".top_Num").css("right", right + 15);
    };
};

$(window).resize(function() {
    var w1 = screen.width;
    var w = document.documentElement.clientWidth;
    var h = document.documentElement.clientHeight;
    var left = (w - minWidth * size) / 2;
    var leftWidth = left - 15;
    var rightWidth = left + 15;
    var shadeLeft = (w - 500) / 2;
    var shadePercent = (w - 500) / 2;
    var shadeTop = (h - 273) / 2;
    var topNumWidth = $(".top_Num").width() + 25;
    if (minWidth * size >= (w - 105)) {
        $(".top_Num").css("right", 77);
    } else if (minWidth * size <= 500) {
        $(".top_Num").css("right", left - topNumWidth);
    } else {
        $(".top_Num").css("right", rightWidth);
    };
    $(".Shadow ").css("width", w - 17);
    $("#body").height(h);
    $(".overflow_ul").css("max-height", h - 100);
    $(".layui-nav").css("max-width", w * 0.29);
    $(".layui-layer").css({
        "left": shadeLeft,
        "top": "60px"
    });
    $(".Tip-page").css({
        "left": shadePercent,
        "top": "60px"
    });
    AddH(size);
});
window.onload = function() {
    parent.parent.postMessage("Page_Loading", "*");
};