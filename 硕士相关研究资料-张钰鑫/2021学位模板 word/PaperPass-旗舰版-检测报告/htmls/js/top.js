(function(System, $) {
    $(function() {

        var treeData = JSON.parse(detail_json);
        tree = loop(treeData);
        $('#detail_report').html(tree);
        $('#detail_report').children('.layui-nav-third-child').removeClass('layui-bg-white').removeClass(
            'layui-nav-third-child');
        var originalData = JSON.parse(original_json);
        original = loopFunc(originalData);
        $("#original_report").html(original);
        $('#original_report').children('.layui-nav-third-child').removeClass('layui-bg-white').removeClass(
            'layui-nav-third-child');
        var proofData = JSON.parse(proof_detail_json);
        proof = loopProof(proofData);
        $('#proof_detail_report').html(proof);
        $('#proof_detail_report').children('.layui-nav-third-child').removeClass('layui-bg-white').removeClass(
            'layui-nav-third-child');

        function loop(dt) {
            var str = "";
            if (dt.length == 0) {
                $(".mlayui-nav-child .mclick_show").hide();
                str += '<div style="color:#000000;text-align:center;margin-left:-40px;padding-top:10px;">' + "目录导航为空" + '</div>';
            } else {
                dt.forEach(item => {
                    var children = (typeof(item.children) != "undefined") && item.children && item.children
                        .length ? loop(item.children) : null;
                    if (children != null) {

                        str +=
                            '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left">' +
                            '<div class="mselect_img select_img_detail"></div>' +
                            '<div class="mselect_tit">' +
                            '<a href="htmls/word/word_report.html' + item.anchorPoint +
                            '" target="left" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' + children +
                            '</div>';
                    } else {

                        str +=
                            '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left ">' +
                            '<div class="mselect_tit leftImg_none">' +
                            '<a href="htmls/word/word_report.html' + item.anchorPoint +
                            '" target="left" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' +
                            '</div>';
                    };
                });
            };
            return str;
        };

        function loopFunc(dt) {
            var str = "";
            if (dt.length == 0) {
                $(".bookmark1_layui-nav-child .mclick_show").hide();
                str += '<div style="color:#000000;text-align:center;margin-left:-40px;padding-top:10px;">' + "目录导航为空" + '</div>';
            } else {
                dt.forEach(item => {
                    var children = (typeof(item.children) != "undefined") && item.children && item.children
                        .length ? loopFunc(item.children) : null;
                    if (children != null) {

                        str += '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left">' +
                            '<div class="mselect_img mselect_img_ori"></div>' +
                            '<div class="mselect_tit">' +
                            '<a href="htmls/word/word_original.html' + item.anchorPoint +
                            '" target="main" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' + children +
                            '</div>';
                    } else {

                        str += '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left ">' +
                            '<div class="mselect_tit leftImg_none">' +
                            '<a href="htmls/word/word_original.html' + item.anchorPoint +
                            '" target="main" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' +
                            '</div>';
                    };
                });
            };
            return str;
        };

        function loopProof(dt) {
            var str = "";
            if (dt.length == 0) {
                $(".proof_layui-nav-child .mclick_show").hide();
                str += '<div style="color:#000000;text-align:center;margin-left:-40px;padding-top:10px;">' + "目录导航为空" + '</div>';
            } else {
                dt.forEach(item => {
                    var children = (typeof(item.children) != "undefined") && item.children && item.children
                        .length ? loopProof(item.children) : null;
                    if (children != null) {

                        str +=
                            '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left">' +
                            '<div class="mselect_img select_img_detail"></div>' +
                            '<div class="mselect_tit">' +
                            '<a href="htmls/word/word_proof.html' + item.anchorPoint +
                            '" target="left" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' + children +
                            '</div>';
                    } else {

                        str +=
                            '<div class="layui-nav-third-child layui-bg-white">' +
                            '<div class="mselect_left ">' +
                            '<div class="mselect_tit leftImg_none">' +
                            '<a href="htmls/word/word_proof.html' + item.anchorPoint +
                            '" target="left" class="third-class">' + item.content + '</a>' +
                            '</div>' +
                            '</div>' +
                            '</div>';
                    };
                });
            };
            return str;
        };


        function PageNav(dt) {
            var str = "";
            if (dt.length == 0) {
                $(".nav-btns").hide();
                str += '<div style="float:none;color:#000000">' + "页面导航为空" + '</div>';
            } else {
                $(".nav-btns").show();
                dt.forEach((element, i) => {
                    str += `<div><a href="javascript:;"><img src="${element}"></a><div class="page-nav-id">${i+1}</div></div>`
                });
            }
            return str;
        }

        function gather(D, num) {
            D.find(".pic-list a").each(function(i, v) {
                var i = i + 1;
                if (i == num) {
                    $(this).parent().parent().find("img").css("border", "1px solid #d6d6d6");
                    D.find(".page-nav-id").css("color", "#343544")
                    $(this).find("img").css("border", "2px solid #35b372");
                    $(this).siblings(".page-nav-id").css("color", "#35b372");
                    $(this).parent()[0].scrollIntoView(false);
                }
            });
            D.find(".pic-list a").click(function() {
                $(this).parent().parent().find("img").css("border", "1px solid #d6d6d6");
                $(this).parent().parent().find(".page-nav-id").css("color", "#343544");
                $(this).find("img").css("border", "2px solid #35b372");
                $(this).siblings(".page-nav-id").css("color", "#35b372");
                var EnterId = $(this).next(".page-nav-id").text();
                var iframe = document.getElementById('childIframe');
                var n_navchild = $(".layui-page-nav-child").css("display");
                var r_navchild = $(".layui-oripage-nav-child").css("display");
                var pt_navchild = $(".layui-pt-nav-child").css("display");
                var opt_navchild = $(".layui-opt-nav-child").css("display");
                var pf_navchild = $(".layui-proof-nav-child").css("display");
                var pft_navchild = $(".layui-proof-pt-nav-child").css("display");
                var obj = {
                    EnterId: EnterId,
                    child3: n_navchild,
                    child4: r_navchild,
                    child5: pt_navchild,
                    child6: opt_navchild,
                    child7: pf_navchild,
                    child8: pft_navchild,
                    page_nav_click: "true"
                }
                iframe.contentWindow.postMessage(obj, "*");
            })
        }

        function scrollPageNum(D, num, click) {
            D.find(".pic-list a").each(function() {
                var EnterId = $(this).next(".page-nav-id").text();
                if (EnterId == num) {
                    $(this).parent().parent().find("img").css("border", "1px solid #d6d6d6");
                    D.find(".page-nav-id").css("color", "#343544")
                    $(this).find("img").css("border", "2px solid #35b372");
                    $(this).siblings(".page-nav-id").css("color", "#35b372");
                    if (click == "false") {
                        $(this).parent()[0].scrollIntoView(false);
                    }
                    var iframe = document.getElementById('childIframe');
                    iframe.contentWindow.postMessage({ page_nav_click: "false" }, "*");
                }
            });
        }

        function FuncPop(tabText) {
            var str = "";
            if (tabText == "详细报告" || tabText == "详细" || tabText == "") {
                var initPage = System.operateState.getInitPage(Report.report_id);
                if (initPage == "switch_text") {
                    str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                } else {
                    str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#484f52 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow: #000000 0px 0px 10px;">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                }
                $(".origina_shade").hide();
                $(".origina_layer").hide();
            } else if (tabText == "文字校审" || tabText == "校审") {
                var proofinitPage = System.operateState.getProofInitPage(Report.report_id);
                if (revise_is_has_word == false || proofinitPage == "proof_switch_text") {
                    str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                } else {
                    str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#484f52 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow: #000000 0px 0px 10px;">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                }
                $(".origina_shade").hide();
                $(".origina_layer").hide();
            } else if (tabText == "report_text" || tabText == "proof_text") {
                str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                $(".main_shade").hide();
                $(".main_layer").hide();
            } else if (tabText == "report_word" || tabText == "proof_word") {
                str += `<div class="layui-layer-shade main_shade" style="z-index: 999; background-color:#484f52 ;display:block;width:70%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt main_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow: #000000 0px 0px 10px;">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载报告，请稍等...</div>
					</div>`;
                $(".main_shade").hide();
                $(".main_layer").hide();
            } else if (tabText == "查看原文" || tabText == "原文" || tabText == "") {
                if (isHasWordReport == true) {
                    var oriState = System.operateState.getOriginalState(Report.report_id);
                    if (oriState == "original_text") {
                        str += `<div class="layui-layer-shade origina_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:100%"></div>
						<div class="layui-layer layui-layer-page layui-layer-prompt origina_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载原文，请稍等...</div>
						</div>`;
                    } else {
                        str += `<div class="layui-layer-shade origina_shade" style="z-index: 999; background-color:#484f52 ;display:block;width:100%"></div>
						<div class="layui-layer layui-layer-page layui-layer-prompt origina_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow: #000000 0px 0px 10px;">
						<div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载原文，请稍等...</div>
						</div>`;
                    }
                } else {
                    str += `<div class="layui-layer-shade origina_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:100%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt origina_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
					   <div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载原文，请稍等...</div>
					</div>`;
                }
                $(".main_shade").hide();
                $(".main_layer").hide();
            } else if (tabText == "view_original") {
                str += `<div class="layui-layer-shade origina_shade" style="z-index: 999; background-color:#f2f2f2 ;display:block;width:100%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt origina_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow:0px 1px 5px 0px rgba(0, 0, 0, 0.11),0px 1px 15px 0px rgba(0, 0, 0, 0.24);">
					   <div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载原文，请稍等...</div>
					</div>`;
                $(".main_shade").hide();
                $(".main_layer").hide();
            } else if (tabText == "word_original") {
                str += `<div class="layui-layer-shade origina_shade" style="z-index: 999; background-color:#484f52 ;display:block;width:100%"></div>
					<div class="layui-layer layui-layer-page layui-layer-prompt origina_layer" style="z-index: 999;display:block;width:400px;height:100px;line-height:100px;box-shadow: #000000 0px 0px 10px;">
					   <div  class="layui-layer-content" style="font-size:18px;padding:0">正在加载原文，请稍等...</div>
					</div>`;
                $(".main_shade").hide();
                $(".main_layer").hide();
            } else {
                $(".main_shade").hide();
                $(".main_layer").hide();
                $(".origina_shade").hide();
                $(".origina_layer").hide();
            };
            $("#m-header").append(str);
            if (tabText == "详细报告" || tabText == "详细" || tabText == "" || tabText == "report_word" || tabText == "report_text"){
                var FoldState = System.operateState.getFoldState(Report.report_id);
                FoldStateWidth(FoldState);
            }
            if(tabText == "文字校审" || tabText == "校审" || tabText == "proof_word" || tabText == "proof_text"){
                var FoldState = System.operateState.getFoldState2(Report.report_id);
                FoldStateWidth(FoldState);
            }
            var mshadeWidth = $(".main_shade").width();
            var shadeLeft = (mshadeWidth - 400) / 2;
            var oshadeWidth = $(".origina_shade").width();
            var originaLeft = (oshadeWidth - 400) / 2;
            $(".main_layer").css("left", shadeLeft);
            $(".origina_layer").css("left", originaLeft);
        };

        function FuncPopRight() {
            var str = '';
            str += `<div class="layui_right_shade">
                        <div class="layui_right_layer">
                            <div class="layui-layer-content" style="font-size:14px;padding:0">正在加载校审结果，请稍等...</div>
                        </div>
                    </div>`;
            $("#m-header").append(str);
        }

        function FoldStateWidth(state){
            if (state == 'none') {
                $('.main_shade').css('width', '100%');
            } else {
                $('.main_shade').css('width', '70%');
            }
        }

        function FunShow() {
            var show = $('.paper-header').css('display');
            if (show == "none") {
                $(".main_shade").css("top", "40px");
                $(".main_layer").css("top", "100px");
                $(".origina_shade").css("top", "40px");
                $(".origina_layer").css("top", "100px");
            } else {
                $(".main_shade").css("top", "140px");
                $(".main_layer").css("top", "200px");
                $(".origina_shade").css("top", "140px");
                $(".origina_layer").css("top", "200px");
            };
        };

        function oriStateFun() {
            var oriState = System.operateState.getOriginalState(Report.report_id);
            if (oriState == "original_text") {
                $(".original-nav a").attr("href", "htmls/view_original.html");
            } else {
                $(".original-nav a").attr("href", "htmls/word/word_original.html");
            }
        }
        var w = document.documentElement.clientWidth;
        navList(w);
        $(".active_li").click(function() {
            var a_Attr = $(this).children("a").attr("href");
            var tab_name = $(this).children("a").text();
            System.operateState.saveALink(Report.report_id, a_Attr);
            System.operateState.saveTabName(Report.report_id, tab_name);
            FuncPop(tab_name);
            FunShow();
            if (tab_name == "查看原文" || tab_name == "原文") {
                if (is_has_word) {
                    oriStateFun();
                }
            }
            if (tab_name == "文字校审" || tab_name == "校审") {
                FuncPopRight();
            } else {
                $(".layui_right_shade").hide();
            }
        });
        var aLink = System.operateState.getALink(Report.report_id);
        var tabName = System.operateState.getTabName(Report.report_id);
        var Typeof = typeof isHasWordReport;
        $("#childIframe").attr("src", aLink);
        $(".active_li").each(function() {
            if ($(this).children("a").text() == tabName) {
                $(this).addClass('active');
                $(this).siblings().removeClass('active');
            };
        });
        if (tabName == "文字校审" || tabName == "校审") {
            FuncPopRight();
        };
        if (Typeof == "boolean" && isHasWordReport == true) {
            FuncPop(tabName);
            oriStateFun();
        };
        if (Typeof == "boolean" && isHasWordReport == false) {
            $(".original-nav a").attr("href", "htmls/view_original.html");
            FuncPop(tabName);
            if (tabName == "详细报告" || tabName == "详细" || tabName == "") {
                FuncPop("report_text");
            }
        };
        

        var html = $(".report_number").html();

        function FuncHeader(Header) {
            var h = document.body.clientHeight;
            var str = "";
            if (Header == 'block' || Header == '') {
                $('.paper-header').css('display', 'none');
                $(".tab").css("padding-right", "0px");
                $(".clickimg").css("background-position", "0 0");
                $("#m-content").height(h - 40);
                $("#childIframe").height(h - 40);
                $(".layui-nav").css("top", "55px");
                $(".overflow_ul").css("height", h - 120);
                $(".main_shade").css("top", "40px");
                $(".main_layer").css("top", "100px");
                $(".origina_shade").css("top", "40px");
                $(".origina_layer").css("top", "100px");
                $(".report_number").css({ "display": "flex", "justify-content": "flex-start", "justify-items": "center" });
                $(".tab_img").attr("src", "htmls/images/Header_logo.png");
                var similarNum = $(".similarNum").html();
                var progress = $(".nav-detail-similar .progress").html();
                var url = $(".report_number a").attr("href");
                str += `<div class="report_number"><a href="${url}" target="_blank" class="g-btn g-btn-default g-btn-sm transparent MLR25" style="margin-top:-4px;">查询真伪</a></div>
					<div style="color:#e6e6e6">总体相似度：</div>
					<div style="width:80px;height:40px;line-height:40px;margin-right:10px;padding-top:4px;">
						<div class="progress" style="margin-top:-5px">${progress}</div>	
					</div>
		            <div style="color:#e6e6e6;display:inline-block">${similarNum}</div>`;
            } else {
                $('.paper-header').css('display', 'block');
                $(".tab").css("padding-right", "113px");
                $(".clickimg").css("background-position", "-56px 0");
                $("#m-content").height(h - 140);
                $("#childIframe").height(h - 140);
                $(".layui-nav").css("top", "155px");
                $(".overflow_ul").css("height", h - 220);
                $(".main_shade").css("top", "140px");
                $(".main_layer").css("top", "200px");
                $(".origina_shade").css("top", "140px");
                $(".origina_layer").css("top", "200px");
                $(".report_number").css("padding-left", "0px");
                $(".tab_img").attr("src", "");
                str = html;
            };
            $(".report_number").html(str);
        };
        var storage = window.localStorage;
        var HeaderState = System.operateState.getHeaderState(Report.report_id);
        FuncHeader(HeaderState);
        $(".HeaderButton").click(function() {
            var show = $('.paper-header').css('display');
            System.operateState.saveHeaderState(Report.report_id, show);
            FuncHeader(show);
        });

        $('.third-class').each(function() {
            var max = 40;
            var str = $(this).text().trim();
            var length = str.length;
            if (length > max) {
                $(this).html(str.substring(0, max) + '...');
            };
        });

        $("#m-nav a").click(function() {
            $(".mlayui-nav-child").hide();
            $(".bookmark1_layui-nav-child").hide();
            $(".proof_layui-nav-child").hide();
            $(".layui-page-nav-child").hide();
            $(".layui-oripage-nav-child").hide();
            $(".layui-pt-nav-child").hide();
            $(".layui-opt-nav-child").hide();
            $(".layui-proof-nav-child").hide();
            $(".layui-proof-pt-nav-child").hide();
            $(".mclick_span").text("全部展开");
            $("i.g-arrow-i").css({
                "border-color": "#56b282 transparent transparent",
                "top": "18px"
            });
        });

        $(".m_close").click(function() {
            $(".mlayui-nav-child").hide();
        });

        $(".mselect_img").on('click', function() {
            if ($(this).hasClass('opened')) {
                $(this).parent().siblings(".layui-nav-third-child").hide();
                $(this).parent().siblings(".layui-nav-third-child").find(".mselect_img").removeClass('opened');
                $(this).parent().siblings(".layui-nav-third-child").find(".mselect_img").css('background-position', '-115px -245px');
                $(this).css('background-position', '-115px -245px');
                $(this).removeClass('opened');
            } else {
                $(this).parent().siblings(".layui-nav-third-child").show();
                $(this).css('background-position', '-77px -245px');
                $(this).addClass('opened');
            };
        });

        $(".select_img_detail").click(function() {
            var flag = true;
            $(".select_img_detail").each(function(k, v) {
                if (!$(v).hasClass('opened')) {
                    flag = false;
                };
            });
            if (flag) {
                $(".mclick_span").text("全部收起");
                $("i.g-arrow-i").css({
                    "border-color": "transparent transparent #56b282",
                    "top": "13px"
                })
            } else {
                $(".mclick_span").text("全部展开");
                $("i.g-arrow-i").css({
                    "border-color": "#56b282 transparent transparent",
                    "top": "18px"
                });
            };
        });

        $(".mselect_img_ori").click(function() {
            var flag = true;
            var mselect_img = $(".mselect_img_ori");
            $.each(mselect_img, function(k, v) {
                if (!$(v).hasClass('opened')) {
                    flag = false;
                };
            });
            if (flag) {
                $(".mclick_span").text("全部收起");
                $("i.g-arrow-i").css({
                    "border-color": "transparent transparent #56b282",
                    "top": "13px"
                });
            } else {
                $(".mclick_span").text("全部展开");
                $("i.g-arrow-i").css({
                    "border-color": "#56b282 transparent transparent",
                    "top": "18px"
                });
            };
        });
        if ($(".layui-nav-third-child").length == 0) {
            $(".mclick_show").hide();
        };

        $(".mclick_show").click(function() {
            var thirdchild = $(this).parent().parent().next().find(".layui-nav-third-child");
            if (thirdchild.is(":hidden")) {
                if ($(".mclick_span").text() == "全部收起") {
                    $(this).parent().parent().next().find(".layui-nav-third-child").hide();
                    $(this).parent().parent().next().find(".mselect_img").css('background-position', '-115px -245px');
                    $(".mclick_span").text("全部展开");
                    $("i.g-arrow-i").css({
                        "border-color": "#56b282 transparent transparent",
                        "top": "18px"
                    });
                    $(".mselect_img").removeClass('opened');
                } else {
                    $(this).parent().parent().next().find(".layui-nav-third-child").show();
                    $(this).parent().parent().next().find(".mselect_img").css('background-position', '-77px -245px');
                    $(".mclick_span").text("全部收起");
                    $("i.g-arrow-i").css({
                        "border-color": "transparent transparent #56b282",
                        "top": "13px"
                    });
                    $(".mselect_img").addClass('opened');
                };
            } else {
                $(this).parent().parent().next().find(".layui-nav-third-child").hide();
                $(this).parent().parent().next().find(".mselect_img").css('background-position', '-115px -245px');
                $(".mclick_span").text("全部展开");
                $("i.g-arrow-i").css({
                    "border-color": "#56b282 transparent transparent",
                    "top": "18px"
                });
                $(".mselect_img").removeClass('opened');
            };
        });

        $(".mselect_tit").mousemove(function() {
            $(this).css("background", "#f2f2f2");
        }).mouseout(function() {
            $(this).css("background", "transparent");
        });

        window.addEventListener("message", function(event) {
            if (event.data == "text") {
                $(".mlayui-nav-child").hide();
                $(".bookmark1_layui-nav-child").hide();
                $(".proof_layui-nav-child").hide();
            };
            if (event.data == "Page_Loading") {
                $(".main_shade").hide();
                $(".main_layer").hide();
                $(".origina_shade").hide();
                $(".origina_layer").hide();
            };
            if (event.data == "Page_Loading_Right") {
                $(".layui_right_shade").hide();
            };
            if (event.data == "Loading_Right") {
                FuncPopRight();
            };
            if (event.data == "report_text") {
                FuncPop("report_text");
                FunShow();
                $(".layui-page-nav-child").hide();
            };
            if (event.data == "report_word") {
                FuncPop("report_word");
                FunShow();
                $(".layui-pt-nav-child").hide();
            };
            if (event.data == "proof_text"){
                FuncPop("proof_text");
                FunShow();
                $(".layui-proof-nav-child").hide();
            };
            if (event.data == "proof_word") {
                FuncPop("proof_word");
                FunShow();
                $(".layui-proof-pt-nav-child").hide();
            };
            if (event.data == "查看原文") {
                FuncPop("查看原文");
                FunShow();
            };
            if (event.data == "view_original") {
                FuncPop("view_original");
                FunShow();
                $(".original-nav a").attr("href", "htmls/view_original.html");
                $(".layui-oripage-nav-child").hide();
                $(".bookmark1_layui-nav-child").hide();
            }
            if (event.data == "word_original") {
                FuncPop("word_original");
                FunShow();
                $(".original-nav a").attr("href", "htmls/word/word_original.html");
                $(".layui-opt-nav-child").hide();
            }
        }, false);

        var iframe = document.getElementById('childIframe');
        var FoldState = System.operateState.getFoldState(Report.report_id);
        var FoldState2 = System.operateState.getFoldState2(Report.report_id);
        window.onload = function() {
            var m_navchild = $(".mlayui-nav-child").css("display");
            var o_navchild = $(".bookmark1_layui-nav-child").css("display");
            var p_navchild = $(".proof_layui-nav-child").css("display");
            var n_navchild = $(".layui-page-nav-child").css("display");
            var r_navchild = $(".layui-oripage-nav-child").css("display");
            var pt_navchild = $(".layui-pt-nav-child").css("display");
            var opt_navchild = $(".layui-opt-nav-child").css("display");
            var pf_navchild = $(".layui-proof-nav-child").css("display");
            var pft_navchild = $(".layui-proof-pt-nav-child").css("display");
            var obj = {
                child0: p_navchild,
                child1: m_navchild,
                child2: o_navchild,
                child3: n_navchild,
                child4: r_navchild,
                child5: pt_navchild,
                child6: opt_navchild,
                child7: pf_navchild,
                child8: pft_navchild,
                FoldState: FoldState,
                FoldState2: FoldState2
            };
            iframe.contentWindow.postMessage(obj, "*");

            $(".m_close").click(function() {
                $(".mlayui-nav-child").hide();
                $(".bookmark1_layui-nav-child").hide();
                $(".proof_layui-nav-child").hide();
                $(".layui-page-nav-child").hide();
                $(".layui-oripage-nav-child").hide();
                $(".layui-pt-nav-child").hide();
                $(".layui-opt-nav-child").hide();
                $(".layui-proof-nav-child").hide();
                $(".layui-proof-pt-nav-child").hide();
                var obj = {
                    child0: p_navchild,
                    child1: m_navchild,
                    child2: o_navchild,
                    child3: n_navchild,
                    child4: r_navchild,
                    child5: pt_navchild,
                    child6: opt_navchild,
                    child7: pf_navchild,
                    child8: pft_navchild,
                };
                iframe.contentWindow.postMessage(obj, "*");
            });
            window.addEventListener("message", function(event) {
                $(".mlayui-nav-child").css("display", event.data.state1);
                $(".bookmark1_layui-nav-child").css("display", event.data.state2);
                $(".proof_layui-nav-child").css("display", event.data.state3);
                $(".layui-page-nav-child").css("display", event.data.navState1);
                $(".layui-oripage-nav-child").css("display", event.data.navState2);
                $(".layui-pt-nav-child").css("display", event.data.navState3);
                $(".layui-opt-nav-child").css("display", event.data.navState4);
                $(".layui-proof-nav-child").css("display", event.data.navState5);
                $(".layui-proof-pt-nav-child").css("display", event.data.navState6);
                if (event.data.FoldState) {
                    BookmarkPosition(event.data.FoldState);
                }
                if (event.data.FoldState2) {
                    BookmarkPosition(event.data.FoldState2);
                }
                if (event.data.section_edit) {
                    $("#section-edit-modal").css("display", "block");
                    $("#layui-layer-section").css("display", "block");
                    var arr = JSON.parse(section_edit_json);
                    var id = event.data.section_edit.trim();
                    $(".section-txt.section-txt-prev").text('');
                    $.each(arr[id], function(i, value) {
                        if (value.score >= 0.4 && value.score < 0.7) {
                            var span = $('<span class="orange">' + value.content + '</span>');
                        } else if (value.score >= 0.7) {
                            var span = $('<span class="red">' + value.content + '</span>');
                        } else {
                            var span = $('<span>' + value.content + '</span>');
                        }
                        $(".section-txt.section-txt-prev").append(span);
                    });
                    $(".section-id").text(id);
                    var textSectionEdit = System.Paper.textSectionEdit();
                    textSectionEdit.call(this);
                }
                if (event.data.modify_document) {
                    $("#modify-document-modal").css("display", "block");
                    $("#layui-layer-document").css("display", "block");
                    $("#modifyDocument").attr("src","htmls/view_original_edit.html");
                }
                if (event.data.page_nav) {
                    pics = PageNav(event.data.page_nav);
                    $("#page_nav .pic-list").html(pics);
                    var DetailNum = System.operateState.getDetailNum(Report.report_id);
                    var num = DetailNum == "" ? "1" : DetailNum;
                    gather($("#page_nav"), num);
                }
                if (event.data.oriPage_nav) {
                    pics = PageNav(event.data.oriPage_nav);
                    $("#oripage_nav .pic-list").html(pics);
                    var OriNum = System.operateState.getOriNum(Report.report_id);
                    var num = OriNum == "" ? "1" : OriNum;
                    gather($("#oripage_nav"), num);
                }
                if (event.data.PageText_nav) {
                    pics = PageNav(event.data.PageText_nav);
                    $("#PageText_nav .pic-list").html(pics);
                    var PlainTextDetailNum = System.operateState.getPlainTextDetailNum(Report.report_id);
                    var num = PlainTextDetailNum == "" ? "1" : PlainTextDetailNum;
                    gather($("#PageText_nav"), num);
                }
                if (event.data.OriPageText_nav) {
                    pics = PageNav(event.data.OriPageText_nav);
                    $("#OriPageText_nav .pic-list").html(pics);
                    var PlainTextOriNum = System.operateState.getPlainTextOriNum(Report.report_id);
                    var num = PlainTextOriNum == "" ? "1" : PlainTextOriNum;
                    gather($("#OriPageText_nav"), num);
                }
                if (event.data.proof_page_nav) {
                    pics = PageNav(event.data.proof_page_nav);
                    $("#proof_page_nav .pic-list").html(pics);
                    var DetailNum = System.operateState.getProofNum(Report.report_id);
                    var num = DetailNum == "" ? "1" : DetailNum;
                    gather($("#proof_page_nav"), num);
                }
                if (event.data.proof_PageText_nav) {
                    pics = PageNav(event.data.proof_PageText_nav);
                    $("#proof_PageText_nav .pic-list").html(pics);
                    var PlainTextDetailNum = System.operateState.getPlainTextProofNum(Report.report_id);
                    var num = PlainTextDetailNum == "" ? "1" : PlainTextDetailNum;
                    gather($("#proof_PageText_nav"), num);
                }

                if (event.data.PageNum1) {
                    var D = $("#page_nav");
                    scrollPageNum(D, event.data.PageNum1, event.data.page_nav_click);
                }
                if (event.data.PageNum2) {
                    var D = $("#oripage_nav");
                    scrollPageNum(D, event.data.PageNum2, event.data.page_nav_click);
                }
                if (event.data.PageNum3) {
                    var D = $("#PageText_nav");
                    scrollPageNum(D, event.data.PageNum3, event.data.page_nav_click);
                }
                if (event.data.PageNum4) {
                    var D = $("#OriPageText_nav");
                    scrollPageNum(D, event.data.PageNum4, event.data.page_nav_click);
                }
                if (event.data.PageNum5) {
                    var D = $("#proof_page_nav");
                    scrollPageNum(D, event.data.PageNum5, event.data.page_nav_click);
                }
                if (event.data.PageNum6) {
                    var D = $("#proof_PageText_nav");
                    scrollPageNum(D, event.data.PageNum6, event.data.page_nav_click);
                }
            }, false);
            var textSectionSave = System.Paper.textSectionSave(null);
            $(document).on('click', '[tpl-section-save=btn]', function() {
                textSectionSave.call(this);
            });
            $(document).on('click', '[tpl-section-cancle=btn],[tpl-section-save=btn],.layui-layer-setwin', function() {
                $("#section-edit-modal").css("display", "none");
                $("#modify-document-modal").css("display", "none");
                $(".layui-layer").css("display", "none");
            });

            var PageNavState = System.operateState.getPageNavState(Report.report_id);
            var OriPageNavState = System.operateState.getOriPageNavState(Report.report_id);
            var PTPageNavState = System.operateState.getPTPageNavState(Report.report_id);
            var OPTPageNavState = System.operateState.getOPTPageNavState(Report.report_id);
            var ProofNavState = System.operateState.getProofNavState(Report.report_id);
            var ProofPtNavState = System.operateState.getProofPtNavState(Report.report_id);
            GetPageNavState($(".layui-page-nav-child"), PageNavState);
            GetPageNavState($(".layui-oripage-nav-child"), OriPageNavState);
            GetPageNavState($(".layui-pt-nav-child"), PTPageNavState);
            GetPageNavState($(".layui-opt-nav-child"), OPTPageNavState);
            GetPageNavState($(".layui-proof-nav-child"), ProofNavState);
            GetPageNavState($(".layui-proof-pt-nav-child"), ProofPtNavState);

            function GetPageNavState(val, num) {
                val.find(".nav-btns li").each(function() {
                    if ($(this).hasClass("active")) {
                        $(this).removeClass("active");
                    }
                });
                if (num == "1") {
                    val.find(".nav-btns li:first-child").addClass("active");
                    val.find(".pic-list").css("width", "144px");
                } else if (num == "2" || num == "") {
                    val.find(".nav-btns li:nth-child(2)").addClass("active");
                    val.find(".pic-list").css("width", "288px");
                } else if (num == "3") {
                    val.find(".nav-btns li:last-child").addClass("active");
                    val.find(".pic-list").css("width", "432px");
                }
            }
            $(document).on('click', '.nav-btns li', function() {
                $(this).addClass("active");
                $(this).siblings().removeClass("active");
                $(this).parents(".layui-nav").find("li").each(function(k, v) {
                    var k = k + 1;
                    if (k == 1 && $(this).hasClass("active")) {
                        $(this).parents(".layui-nav").find(".pic-list").css("width", "144px");
                        SaveState($(this), k);
                    } else if (k == 2 && $(this).hasClass("active")) {
                        $(this).parents(".layui-nav").find(".pic-list").css("width", "288px");
                        SaveState($(this), k);
                    } else if (k == 3 && $(this).hasClass("active")) {
                        $(this).parents(".layui-nav").find(".pic-list").css("width", "432px");
                        SaveState($(this), k);
                    }
                });
            });

            function SaveState(d, k) {
                var pageNav = d.parents(".layui-bg-gray").hasClass("layui-page-nav-child");
                var oriPageNav = d.parents(".layui-bg-gray").hasClass("layui-oripage-nav-child");
                var ptPageNav = d.parents(".layui-bg-gray").hasClass("layui-pt-nav-child");
                var optPageNav = d.parents(".layui-bg-gray").hasClass("layui-opt-nav-child");
                var proofNav = d.parents(".layui-bg-gray").hasClass("layui-proof-nav-child");
                var proofPtNav = d.parents(".layui-bg-gray").hasClass("layui-proof-pt-nav-child");
                if (pageNav) {
                    System.operateState.savePageNavState(Report.report_id, k);
                } else if (oriPageNav) {
                    System.operateState.saveOriPageNavState(Report.report_id, k);
                } else if (ptPageNav) {
                    System.operateState.savePTPageNavState(Report.report_id, k);
                } else if (optPageNav) {
                    System.operateState.saveOPTPageNavState(Report.report_id, k);
                } else if (proofNav) {
                    System.operateState.saveProofNavState(Report.report_id, k);
                } else if (proofPtNav) {
                    System.operateState.saveProofPtNavState(Report.report_id, k);
                }
            }
        };

        var w1 = screen.width;
        var w = document.documentElement.clientWidth;
        var h = document.documentElement.clientHeight;
        var tabName = System.operateState.getTabName(Report.report_id);
        var mshadeWidth = $(".main_shade").width();
        var shadeLeft = (mshadeWidth - 400) / 2;
        var oshadeWidth = $(".origina_shade").width();
        var originaLeft = (oshadeWidth - 400) / 2;
        $("body").height(h);
        $(".main_layer").css("left", shadeLeft);
        $(".origina_layer").css("left", originaLeft);
        if(tabName == "" || tabName == "详细报告" || tabName == "详细") {
            var FoldState = System.operateState.getFoldState(Report.report_id);
            BookmarkPosition(FoldState);
        }else if(tabName == "文字校审" || tabName == "校审"){
            var FoldState2 = System.operateState.getFoldState2(Report.report_id);
            BookmarkPosition(FoldState2);
        }
        $(".bookmark1").css({
            "width": w * 0.29,
            "right": 72
        });
        listHeight();

        function BookmarkPosition(FoldState) {
            var w = document.documentElement.clientWidth;
            if (FoldState == 'none') {
                $(".bookmark").css({
                    "width": w * 0.29,
                    "left": "auto",
                    "right": "72px"
                });
            } else {
                $(".bookmark").css({
                    "width": w * 0.29,
                    "left": w * 0.7 - 30,
                    "right": "auto"
                });
            }
        }
    });

    function listHeight() {
        var h = document.documentElement.clientHeight;
        var show = $('.paper-header').css('display');
        if (show == 'block') {
            $(".overflow_ul").css("height", h - 220);
        };
        if (show == 'none') {
            $(".overflow_ul").css("height", h - 120);
        };
    };

    function navList(w) { 
        $("#nav_ul li").each(function(i,v){
            var text = ["详细报告","文字校审","综合评估","查看原文","打印PDF","使用帮助"];
            var txt = ["详细","校审","评估","原文","PDF","帮助"];
            var curItem = $(this).find("a").html();
            if (w > 1280) {
                $(this).find("a").html(text[txt.indexOf(curItem)]);
                $(".nav_tip").remove();
                $('.tab li a').css("padding", "0 20px");
                $(".nav-tit h4").css("font-size", "20px");
            } else {
                $(this).find("a").html(txt[text.indexOf(curItem)]);
                if ($(this).find(".nav_tip").length < 1) { 
                    $(this).append('<div class="nav_tip">'+ text[text.indexOf(curItem)] +'</div>');
                }
                $('.tab li a').css("padding", "0 10px");
                $(".nav-tit h4").css("font-size", "16px");
            }
        })
        $("#nav_ul li").mousemove(function() {
            $(this).find(".nav_tip").css("display", "block");
        }).mouseout(function() {
            $(this).find('.nav_tip').css("display", "none");
        });
    }

    $(window).resize(function() {
        var w1 = screen.width;
        var w = document.documentElement.clientWidth;
        var h = document.documentElement.clientHeight;
        var tabName = System.operateState.getTabName(Report.report_id);
        var mshadeWidth = $(".main_shade").width();
        var shadeLeft = (mshadeWidth - 400) / 2;
        var oshadeWidth = $(".origina_shade").width();
        var originaLeft = (oshadeWidth - 400) / 2;
        $("body").height(h);
        $(".main_layer").css("left", shadeLeft);
        $(".origina_layer").css("left", originaLeft);
        var FoldState = System.operateState.getFoldState(Report.report_id);
        if(tabName == "" || tabName == "详细报告" || tabName == "详细") {
            var FoldState = System.operateState.getFoldState(Report.report_id);
            BookmarkPosition(FoldState);
        }else if(tabName == "文字校审" || tabName == "校审"){
            var FoldState2 = System.operateState.getFoldState2(Report.report_id);
            BookmarkPosition(FoldState2);
        }
        function BookmarkPosition(FoldState) {
            var w = document.documentElement.clientWidth;
            if (FoldState == 'block' || FoldState == '') {
                $(".bookmark").css({
                    "width": w * 0.29,
                    "left": w * 0.7 - 30,
                    "right": "auto"
                });
            } else {
                $(".bookmark").css({
                    "width": w * 0.29,
                    "left": "auto",
                    "right": "72px"
                });
            }
        }
        $(".bookmark1").css({
            "width": w * 0.29,

            "right": 72
        });
        listHeight();
        navList(w);
    });
})(Report, jQuery);