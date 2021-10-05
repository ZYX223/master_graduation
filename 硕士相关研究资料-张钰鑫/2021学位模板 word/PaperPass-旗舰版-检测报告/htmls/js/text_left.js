(function(System, $) {
    $(function() {



        parent.parent.postMessage("word", "*");
        var storage = window.localStorage;
        var LineState = System.operateState.getPlainTextLineState(Report.report_id);
        var TypeOf = typeof isHasDetailReport;
        var ProTypeOf = typeof isHasProofReport;
        if (TypeOf == "boolean" && isHasDetailReport == true) {
            var DetailNum = System.operateState.getPlainTextDetailNum(Report.report_id);
            var Detailoffset = System.operateState.getPlainTextDetailoffset(Report.report_id);
            var DetailSize = System.operateState.getPlainTextDetailSize(Report.report_id);
            var DetailTopNum = System.operateState.getPlainTextDetailTopNum(Report.report_id);
            var DetailPercent = System.operateState.getPlainTextDetailPercent(Report.report_id);
            var ZoomBtnState = System.operateState.getPlainTextZoomBtnState(Report.report_id);
            var FoldState = System.operateState.getFoldState(Report.report_id);
        } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
            var ProofNum = System.operateState.getPlainTextProofNum(Report.report_id);
            var ProofOffset = System.operateState.getPlainTextProofOffset(Report.report_id);
            var ProofSize = System.operateState.getPlainTextProofSize(Report.report_id);
            var ProofTopNum = System.operateState.getPlainTextProofTopNum(Report.report_id);
            var ProofPercent = System.operateState.getPlainTextProofPercent(Report.report_id);
            var ProZoomBtnState = System.operateState.getPlainTextProZoomBtnState(Report.report_id);
            var FoldState2 = System.operateState.getFoldState2(Report.report_id);
        } else {
            var OriNum = System.operateState.getPlainTextOriNum(Report.report_id);
            var Orioffset = System.operateState.getPlainTextOrioffset(Report.report_id);
            var OriSize = System.operateState.getPlainTextOriSize(Report.report_id);
            var OriTopNum = System.operateState.getPlainTextOriTopNum(Report.report_id);
            var OrilPercent = System.operateState.getPlainTextOrilPercent(Report.report_id);
            var OriZoomBtnState = System.operateState.getPlainTextOriZoomBtnState(Report.report_id);
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
                $("#plain-text-page").scrollTop(offsetH);
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
                $("#menuButton").css("background", "#f2f2f2");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();

            };




        } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
            parentWidth(ProofSize ? ProofSize : 1);
            AddH(ProofSize ? ProofSize : 1);
            FunPage(ProofNum, ProofOffset);
            size = (ProofSize ? ProofSize : 1);
            if (ProZoomBtnState == "" || ProZoomBtnState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#f2f2f2");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
            };
        } else {
            parentWidth(OriSize ? OriSize : 1);
            AddH(OriSize ? OriSize : 1);
            FunPage(OriNum, Orioffset);
            size = (OriSize ? OriSize : 1);

            if (OriZoomBtnState == "" || OriZoomBtnState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#f2f2f2");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();

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
            $(".layui-layer-content").html("确认切换至“Word模式”吗？");

            if (TypeOf == "boolean" && isHasDetailReport == true) {
                $(".detail_btn").data('type', 'switch_word');
            }
            if (ProTypeOf == "boolean" && isHasProofReport == true) {
                $(".detail_btn").data('type', 'proof_switch_word');
            }
        });
        $(".detail_btn").click(function() {
            type = $(this).data('type');
            if (type == 'switch_word') {
                $(".layui-layer-btn0").attr("href", $(this).data('txtLink'));
                parent.parent.postMessage("report_word", "*");
                parent.postMessage("text_right", "*");
                System.operateState.saveInitPage(Report.report_id, type);
            }
            if (type == 'proof_switch_word') {
                $(".layui-layer-btn0").attr("href", $(this).data('txtLink'));
                parent.postMessage("proofPt_right", "*");
                parent.parent.postMessage("proof_word", "*");
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
				$(".subline").css("background-position","-161px -9px");
			} else {
				$(".red").parent().addClass('red_a');
				$(".orange").parent().addClass('orange_a');
				$(".subline").css("background-position","-43px -153px");
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
				System.operateState.savePlainTextLineState(Report.report_id, str_has);
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
                        PageNum3: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.savePlainTextDetailNum(Report.report_id, Page);
                    System.operateState.savePlainTextDetailoffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.savePlainTextDetailSize(Report.report_id, PageSize);
                System.operateState.savePlainTextDetailPercent(Report.report_id, percent);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                $("body").off("scroll").on("scroll", function(e) {
                    scrollFunc(_Size);
                    var offsetH = $("body").scrollTop();
                    var Page = nowPage;
                    var obj = {
                        PageNum6: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.savePlainTextProofNum(Report.report_id, Page);
                    System.operateState.savePlainTextProofOffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.savePlainTextProofSize(Report.report_id, PageSize);
                System.operateState.savePlainTextProofPercent(Report.report_id, percent);
            } else {
                $("body").off("scroll").on("scroll", function(e) {
                    scrollFunc(_Size);
                    var offsetH = $("body").scrollTop();
                    var Page = nowPage;
                    var obj = {
                        PageNum4: Page,
                        page_nav_click: page_nav_click
                    }
                    parent.parent.postMessage(obj, "*");
                    System.operateState.savePlainTextOriNum(Report.report_id, Page);
                    System.operateState.savePlainTextOrioffset(Report.report_id, offsetH);
                });
                var PageSize = _Size;
                System.operateState.savePlainTextOriSize(Report.report_id, PageSize);
                System.operateState.savePlainTextOrilPercent(Report.report_id, percent);
            };
        };

        function TopNum() {
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                var TopNum = $(".top_Num").css("right");
                System.operateState.savePlainTextDetailTopNum(Report.report_id, TopNum);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                var TopNum = $(".top_Num").css("right");
                System.operateState.savePlainTextProofTopNum(Report.report_id, TopNum);
            } else {
                var TopNum = $(".top_Num").css("right");
                System.operateState.savePlainTextOriTopNum(Report.report_id, TopNum);
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
            $(this).css("background-image", "-webkit-linear-gradient( 90deg, rgb(216,216,216) 0%, rgb(239,239,239) 100%)");
            $(this).children(".bot_tip").css("display", "block");
        }).mouseout(function() {
            $(this).css("background", "#f2f2f2");
            $(this).children(".bot_tip").css("display", "none");
        });

        $("#CollagenButton").click(function(e) {
            var listState = $(".percent_list").css("display");
            if (listState == "block") {
                $(".percent_list").hide();
                $(this).children(".percent_tip").css("display", "block");
            } else {
                $(".percent_list").show();
                $("#CollagenButton").css("background-image", "-webkit-linear-gradient( 90deg, rgb(224,224,224) 0%, rgb(170,170,170) 170%)");
            };
            $(this).children(".percent_tip").css("display", "none");
            e.stopPropagation();
        }).mousemove(function() {
            var listState = $(".percent_list").css("display");
            if (listState == "block") {
                $(this).children(".percent_tip").css("display", "none");
                $("#CollagenButton").css("background-image", "-webkit-linear-gradient( 90deg, rgb(224,224,224) 0%, rgb(170,170,170) 170%)");
            }
        }).mouseout(function() {
            if ($(".percent_list").css("display") == "block") {
                $("#CollagenButton").css("background-image", "-webkit-linear-gradient( 90deg, rgb(224,224,224) 0%, rgb(170,170,170) 170%)");
            } else {
                $("#CollagenButton").css("background", "#f2f2f2");
            };
            $(this).children(".percent_tip").css("display", "none");
        });
        $(document).click(function() {
            $(".percent_list").hide();
            $("#CollagenButton").css("background", "#f2f2f2");
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
            $(this).css("background", "#e9e9e9");
        }).mouseout(function() {
            $(this).css("background", "#f2f2f2");
        });;

        $("#menuButton").click(function() {


            var showState = $(".allButton").css("display");
            $(".percent_list").hide();
            $("#CollagenButton").css("background", "#f2f2f2");
            if (TypeOf == "boolean" && isHasDetailReport == true) {
                var showState = $(".allButton").css("display");
                System.operateState.savePlainTextZoomBtnState(Report.report_id, showState);
            } else if (ProTypeOf == "boolean" && isHasProofReport == true) {
                var showState = $(".allButton").css("display");
                System.operateState.savePlainTextProZoomBtnState(Report.report_id, showState);
            } else {
                var showState = $(".allButton").css("display");
                System.operateState.savePlainTextOriZoomBtnState(Report.report_id, showState);
            };
            if (showState == "block") {
                $(".allButton").hide();
                $(".CollagenButton").hide();
                $("#menuButton").css("background", "#f2f2f2");
            } else {
                $(".allButton").show();
                $(".CollagenButton").show();
            };








        }).mousemove(function() {
            $(this).children(".menu_tip").css("display", "block");
        }).mouseout(function() {
            $(this).children(".menu_tip").css("display", "none");
            $("#menuButton").css("background", "#f2f2f2");
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
            $(".input_pageNum").css("background-color", "#fbfbfb");
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
            $(".input_pageNum").css("background-color", "#fbfbfb");
        }).mouseout(function() {
            if ($('.input_pageNum').is(':focus')) {
                $(".input_pageNum").css("background-color", "#fbfbfb");
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
                "border": "1px solid #dbdbdb",
                "background-color": "#fbfbfb"
            });
            $(this).children(".tip").css("display", "block");
        }).mouseout(function() {
            $(this).css({
                "background-color": "transparent",
                "border": "1px solid transparent"
            });
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
        $("#plain-text-page").height(h);
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
                    var _size = System.operateState.getPlainTextDetailSize(Report.report_id);
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

        $(document).on('click', '.edit_icon', function() {
            var id = $(this).prev().text();
            var obj = { "section_edit": id };
            parent.parent.postMessage(obj, "*");
        })

        $(document).on('click', '.modify_document', function() {
            var obj = { "modify_document": true }
            parent.parent.postMessage(obj, "*");
        })
        $(".ori_switch_btn").click(function() {
            parent.parent.postMessage("word_original", "*");
            System.operateState.saveOriginalState(Report.report_id, "original_word");
            System.operateState.saveALink(Report.report_id, "htmls/word/word_original.html");
        })
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
    $("#plain-text-page").height(h);
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