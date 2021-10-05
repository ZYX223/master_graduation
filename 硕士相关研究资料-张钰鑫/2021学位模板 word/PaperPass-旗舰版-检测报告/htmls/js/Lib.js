

if (!GRN_P) { var GRN_P = 'Report'; }
(function(global, namespace, factory, undefined) {
    'use strict';
    global[namespace] = factory(global, namespace);
})(typeof window !== "undefined" ? window : this, GRN_P, function(W, namespace, undefined) {
    'use strict';
    return {};
});

(function(System, $) {
    var __this__ = null;

    function isset(s) {
        return (typeof s !== "undefined" && s !== null);
    }
    var cache = [],
        cache_name = 'cache',
        myStorage = null;

    function isStorage() {
        if (typeof(myStorage) !== "undefined") {
            return true;
        } else {
            return false;
        }
    }

    function set(value, name) {
        value = value || cache;
        name = name || cache_name;
        if (isStorage()) {
            remove();
            myStorage.setItem(name, JSON.stringify(value));
        }
    }

    function remove(name) {
        name = name || cache_name;
        if (isStorage()) {
            myStorage.removeItem(name);
        }
    }

    function clear() {
        if (isStorage()) {
            myStorage.clear();
        } else {
            cache = [];
        }

    }

    function get(name) {
        name = name || cache_name;
        if (isStorage()) {
            return (JSON.parse(myStorage.getItem(name))) || cache;
        } else {
            return cache;
        }

    }

    
    function Cache(name, type) {
        __this__ = this;
        this.caches = [];
        cache_name = name || 'cache';
        this.cache_name = cache_name;
        this.myStorage = type || localStorage;
    }
    Cache.prototype = {
        'constructor': Cache,
        '_className': 'Cache',
        'init': function() {
            myStorage = this.myStorage;
            cache_name = this.cache_name;
            cache = get();
        },
        'cache': function(key, value, callback) {
            this.init();
            cache = get();
            var index = this.exist(key, value);
            if ($.isFunction(callback)) {
                callback.call(this, index);
            }
            return index;
        },
        'set': function(Obj, key, value) {
            this.init();
            Obj[key] = value;
            cache.push(Obj);
            set();

        },
        'update': function(index, Obj) {
            this.init();
            cache[index] = Obj;
            set();
        },
        'get': function(index) {
            this.init();
            return isset(index) ? cache[index] : cache;
        },
        'exist': function(key, value) {
            this.init();
            for (var i = 0, len = cache.length; i < len; i++) {
                if ((key in cache[i]) && (value == cache[i][key])) {
                    return i;
                }
            }
            return -1;
        },

        'clear': function() {
            this.init();
            this.remove();
            clear();
        },
        'remove': function(index) {
            this.init();
            if (index) {
                if (index > -1 && index < cache.length - 1) {
                    cache.splice(index, 1);

                }
            } else {
                cache = [];
            }
            set();
        }
    };
    System.Cache = Cache;
})(Report, jQuery);



(function(System, $) {
    'use strict';
    var storage = window.localStorage;

    function getObj(id) {
        id = id + "_s";
        var obj = JSON.parse(storage.getItem(id))
        return obj
    }

    function setObj(id, obj) {
        id = id + "_s";
        storage.setItem(id, JSON.stringify(obj))
    }

    function operateState() {}
    
    operateState.saveId = function(id) {
        var defaults = {
            "red": "",
            "LineState": "",
            "HeaderState": '',
            "aLink": "htmls/detail_report.html",
            "tabName": "",
            "DetailNum": "",
            "Detailoffset": "",
            "DetailSize": "",
            "DetailTopNum": "",
            'DetailPercent': "",
            "ZoomBtnState": "",
            "OriNum": "",
            "Orioffset": "",
            "OriSize": "",
            "OriTopNum": "",
            'OrilPercent': "",
            "OriZoomBtnState": "",
            "FoldState": "",
            "InitPage": "",
            "PlainTextLineState": "",
            "PlainTextDetailNum": "",
            "PlainTextDetailoffset": "",
            "PlainTextDetailSize": "",
            "PlainTextDetailPercent": "",
            "PlainTextDetailTopNum": "",
            "PlainTextZoomBtnState": "",
            "PlainTextOriNum": "",
            "PlainTextOrioffset": "",
            "PlainTextOriSize": "",
            "PlainTextOriTopNum": "",
            'PlainTextOrilPercent': "",
            "PlainTextOriZoomBtnState": "",
            "OriginalState": "",
            "PageNavState": "",
            "OriPageNavState": "",
            "PTPageNavState": "",
            "OPTPageNavState": "",

            "ProofNum": "",
            "ProofOffset": "",
            "ProofSize": "",
            "ProofTopNum": "",
            "ProofPercent": "",
            "ProZoomBtnState": "",
            "ProofNavState": "",
            "ProofInitPage": "",
            "FoldState2": "",
            "PlainTextProofNum": "",
            "PlainTextProofOffset": "",
            "PlainTextProofSize": "",
            "PlainTextProofTopNum": "",
            "PlainTextProofPercent": "",
            "PlainTextProZoomBtnState": "",
            "ProofPtNavState": "",
            "Suggestions": [],
            "ProofTab": "",
        };
        id = id + "_s";
        if (storage.getItem(id) === null) {
            storage.setItem(id, JSON.stringify(defaults));
        } else {
            return
        }

    };
    
    operateState.saveHeaderState = function(id, show) {
        var obj = getObj(id);
        obj.HeaderState = show;
        setObj(id, obj);
    };
    operateState.getHeaderState = function(id) {
        var obj = getObj(id);
        return obj.HeaderState;
    };
    
    operateState.saveALink = function(id, a_Attr) {
        var obj = getObj(id);
        obj.aLink = a_Attr;
        setObj(id, obj);
    };
    operateState.getALink = function(id) {
        var obj = getObj(id);
        return obj.aLink;
    };
    
    operateState.saveTabName = function(id, tab_name) {
        var obj = getObj(id);
        obj.tabName = tab_name;
        setObj(id, obj);
    };
    operateState.getTabName = function(id) {
        var obj = getObj(id);
        return obj.tabName;
    };
    
    operateState.saveLineState = function(id, isUsed) {
        var obj = getObj(id);
        obj.LineState = isUsed;
        setObj(id, obj);
    };
    operateState.getLineState = function(id) {
        var obj = getObj(id);
        return obj.LineState;
    };
    
    operateState.saveDetailNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.DetailNum = PageNum;
        setObj(id, obj);
    };
    operateState.getDetailNum = function(id) {
        var obj = getObj(id);
        return obj.DetailNum;
    };
    
    operateState.saveDetailTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.DetailTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getDetailTopNum = function(id) {
        var obj = getObj(id);
        return obj.DetailTopNum;
    };
    
    operateState.saveDetailoffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.Detailoffset = offsetH;
        setObj(id, obj);
    };
    operateState.getDetailoffset = function(id) {
        var obj = getObj(id);
        return obj.Detailoffset;
    };
    
    operateState.saveDetailSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.DetailSize = PageSize;
        setObj(id, obj);
    };
    operateState.getDetailSize = function(id) {
        var obj = getObj(id);
        return obj.DetailSize;
    };
    
    operateState.saveDetailPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.DetailPercent = Percent;
        setObj(id, obj);
    };
    operateState.getDetailPercent = function(id) {
        var obj = getObj(id);
        return obj.DetailPercent;
    };
    
    operateState.saveZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.ZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.ZoomBtnState;
    };

    
    operateState.saveOriNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.OriNum = PageNum;
        setObj(id, obj);
    };
    operateState.getOriNum = function(id) {
        var obj = getObj(id);
        return obj.OriNum;
    };
    
    operateState.saveOrioffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.Orioffset = offsetH;
        setObj(id, obj);

    };
    operateState.getOrioffset = function(id) {
        var obj = getObj(id);
        return obj.Orioffset;
    };
    
    operateState.saveOriTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.OriTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getOriTopNum = function(id) {
        var obj = getObj(id);
        return obj.OriTopNum;
    };

    
    operateState.saveOriSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.OriSize = PageSize;
        setObj(id, obj);
    };
    operateState.getOriSize = function(id) {
        var obj = getObj(id);
        return obj.OriSize;
    };
    
    operateState.saveOrilPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.OrilPercent = Percent;
        setObj(id, obj);
    };
    operateState.getOrilPercent = function(id) {
        var obj = getObj(id);
        return obj.OrilPercent;
    };
    
    operateState.saveOriZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.OriZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getOriZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.OriZoomBtnState;
    };
    
    operateState.saveFoldState = function(id, isFold) {
        var obj = getObj(id);
        obj.FoldState = isFold;
        setObj(id, obj);
    };
    operateState.getFoldState = function(id) {
        var obj = getObj(id);
        return obj.FoldState;
    };

    
    operateState.saveInitPage = function(id, page) {
        var obj = getObj(id);
        obj.InitPage = page;
        setObj(id, obj);
    };
    operateState.getInitPage = function(id) {
        var obj = getObj(id);
        return obj.InitPage;
    };
    
    operateState.savePlainTextLineState = function(id, isUsed) {
        var obj = getObj(id);
        obj.PlainTextLineState = isUsed;
        setObj(id, obj);
    };
    operateState.getPlainTextLineState = function(id) {
        var obj = getObj(id);
        return obj.PlainTextLineState;
    };
    
    operateState.savePlainTextDetailNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.PlainTextDetailNum = PageNum;
        setObj(id, obj);
    };
    operateState.getPlainTextDetailNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextDetailNum;
    };
    
    operateState.savePlainTextDetailoffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.PlainTextDetailoffset = offsetH;
        setObj(id, obj);
    };
    operateState.getPlainTextDetailoffset = function(id) {
        var obj = getObj(id);
        return obj.PlainTextDetailoffset;
    };
    
    operateState.savePlainTextDetailSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.PlainTextDetailSize = PageSize;
        setObj(id, obj);
    };
    operateState.getPlainTextDetailSize = function(id) {
        var obj = getObj(id);
        return obj.PlainTextDetailSize;
    };
    
    operateState.savePlainTextDetailPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.PlainTextDetailPercent = Percent;
        setObj(id, obj);
    };
    operateState.getPlainTextDetailPercent = function(id) {
        var obj = getObj(id);
        return obj.PlainTextDetailPercent;
    };
    
    operateState.savePlainTextDetailTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.PlainTextDetailTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getPlainTextDetailTopNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextDetailTopNum;
    };
    
    operateState.savePlainTextZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.PlainTextZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getPlainTextZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.PlainTextZoomBtnState;
    };
    
    operateState.savePlainTextOriNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.PlainTextOriNum = PageNum;
        setObj(id, obj);
    };
    operateState.getPlainTextOriNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOriNum;
    };
    
    operateState.savePlainTextOrioffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.PlainTextOrioffset = offsetH;
        setObj(id, obj);

    };
    operateState.getPlainTextOrioffset = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOrioffset;
    };
    
    operateState.savePlainTextOriTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.PlainTextOriTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getPlainTextOriTopNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOriTopNum;
    };
    
    operateState.savePlainTextOriSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.PlainTextOriSize = PageSize;
        setObj(id, obj);
    };
    operateState.getPlainTextOriSize = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOriSize;
    };
    
    operateState.savePlainTextOrilPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.PlainTextOrilPercent = Percent;
        setObj(id, obj);
    };
    operateState.getPlainTextOrilPercent = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOrilPercent;
    };
    
    operateState.savePlainTextOriZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.PlainTextOriZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getPlainTextOriZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.PlainTextOriZoomBtnState;
    };
    
    operateState.saveOriginalState = function(id, isShow) {
        var obj = getObj(id);
        obj.OriginalState = isShow;
        setObj(id, obj);
    };
    operateState.getOriginalState = function(id) {
        var obj = getObj(id);
        return obj.OriginalState;
    };
    
    operateState.savePageNavState = function(id, num) {
        var obj = getObj(id);
        obj.PageNavState = num;
        setObj(id, obj);
    }
    operateState.getPageNavState = function(id) {
            var obj = getObj(id);
            return obj.PageNavState;
    }
    
    operateState.saveOriPageNavState = function(id, num) {
        var obj = getObj(id);
        obj.OriPageNavState = num;
        setObj(id, obj);
    }
    operateState.getOriPageNavState = function(id) {
            var obj = getObj(id);
            return obj.OriPageNavState;
        }
    
    operateState.savePTPageNavState = function(id, num) {
        var obj = getObj(id);
        obj.PTPageNavState = num;
        setObj(id, obj);
    }
    operateState.getPTPageNavState = function(id) {
            var obj = getObj(id);
            return obj.PTPageNavState;
        }
    
    operateState.saveOPTPageNavState = function(id, num) {
        var obj = getObj(id);
        obj.OPTPageNavState = num;
        setObj(id, obj);
    }
    operateState.getOPTPageNavState = function(id) {
        var obj = getObj(id);
        return obj.OPTPageNavState;
    }

    
    operateState.saveProofNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.ProofNum = PageNum;
        setObj(id, obj);
    };
    operateState.getProofNum = function(id) {
        var obj = getObj(id);
        return obj.ProofNum;
    };
    
    operateState.saveProofTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.ProofTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getProofTopNum = function(id) {
        var obj = getObj(id);
        return obj.ProofTopNum;
    };
    
    operateState.saveProofOffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.ProofOffset = offsetH;
        setObj(id, obj);
    };
    operateState.getProofOffset = function(id) {
        var obj = getObj(id);
        return obj.ProofOffset;
    };
    
    operateState.saveProofSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.ProofSize = PageSize;
        setObj(id, obj);
    };
    operateState.getProofSize = function(id) {
        var obj = getObj(id);
        return obj.ProofSize;
    };
    
    operateState.saveProofPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.ProofPercent = Percent;
        setObj(id, obj);
    };
    operateState.getProofPercent = function(id) {
        var obj = getObj(id);
        return obj.ProofPercent;
    };
    
    operateState.saveProZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.ProZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getProZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.ProZoomBtnState;
    };

    
    operateState.savePlainTextProofNum = function(id, PageNum) {
        var obj = getObj(id);
        obj.PlainTextProofNum = PageNum;
        setObj(id, obj);
    };
    operateState.getPlainTextProofNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProofNum;
    };
    
    operateState.savePlainTextProofTopNum = function(id, TopNum) {
        var obj = getObj(id);
        obj.PlainTextProofTopNum = TopNum;
        setObj(id, obj);
    };
    operateState.getPlainTextProofTopNum = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProofTopNum;
    };
    
    operateState.savePlainTextProofOffset = function(id, offsetH) {
        var obj = getObj(id);
        obj.PlainTextProofOffset = offsetH;
        setObj(id, obj);
    };
    operateState.getPlainTextProofOffset = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProofOffset;
    };
    
    operateState.savePlainTextProofSize = function(id, PageSize) {
        var obj = getObj(id);
        obj.PlainTextProofSize = PageSize;
        setObj(id, obj);
    };
    operateState.getPlainTextProofSize = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProofSize;
    };
    
    operateState.savePlainTextProofPercent = function(id, Percent) {
        var obj = getObj(id);
        obj.PlainTextProofPercent = Percent;
        setObj(id, obj);
    };
    operateState.getPlainTextProofPercent = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProofPercent;
    };
    
    operateState.savePlainTextProZoomBtnState = function(id, isShow) {
        var obj = getObj(id);
        obj.PlainTextProZoomBtnState = isShow;
        setObj(id, obj);
    };
    operateState.getPlainTextProZoomBtnState = function(id) {
        var obj = getObj(id);
        return obj.PlainTextProZoomBtnState;
    };
    
    operateState.saveProofNavState = function(id, num) {
        var obj = getObj(id);
        obj.ProofNavState = num;
        setObj(id, obj);
    }
    operateState.getProofNavState = function(id) {
            var obj = getObj(id);
            return obj.ProofNavState;
    }
    
    operateState.saveProofPtNavState = function(id, num) {
        var obj = getObj(id);
        obj.ProofPtNavState = num;
        setObj(id, obj);
    }
    operateState.getProofPtNavState = function(id) {
            var obj = getObj(id);
            return obj.ProofPtNavState;
    }
    
    operateState.saveProofInitPage = function(id, page) {
        var obj = getObj(id);
        obj.ProofInitPage = page;
        setObj(id, obj);
    };
    operateState.getProofInitPage = function(id) {
        var obj = getObj(id);
        return obj.ProofInitPage;
    };
    
    operateState.saveSuggestions = function(id, suggestions) {
        var obj = getObj(id);
        obj.Suggestions = suggestions;
        setObj(id, obj);
    };
    operateState.getSuggestions = function(id) {
        var obj = getObj(id);
        return obj.Suggestions;
    };
    
    operateState.saveProofTab = function(id, type) {
        var obj = getObj(id);
        obj.ProofTab = type;
        setObj(id, obj);
    };
    operateState.getProofTab = function(id) {
        var obj = getObj(id);
        return obj.ProofTab;
    };
    
    operateState.saveFoldState2 = function(id, isFold) {
        var obj = getObj(id);
        obj.FoldState2 = isFold;
        setObj(id, obj);
    };
    operateState.getFoldState2 = function(id) {
        var obj = getObj(id);
        return obj.FoldState2;
    };

    System.operateState = operateState
})(Report, jQuery);



(function(System, $) {
    'use strict';
    var WIN_H, MIN_HEIGHT = "600";
    var cache = null;

    function getCache() {
        if (!cache) { cache = new System.Cache(System.report_id, localStorage); }
    }

    function Paper() {}
    
    Paper.setMainIframeHeight = function(D) {
        var defaults = {
            "header": '#m-header',
            "nav": '#m-nav',
            "content": '#m-content'
        };
        D = $.isPlainObject(D) ? $.extend(defaults, D) : defaults;
        var $iframe, $header, header_h, $nav, nav_h;

        function initIframe() {
            var height = WIN_H - (header_h + nav_h);




            $iframe.height(height);
            $iframe.find('iframe').height(height);
        }
        return function() {
            $header = $(D.header);
            header_h = $header.height();
            if (!header_h) { header_h = 0; }
            $nav = $(D.nav);
            nav_h = $nav.height();
            if (!nav_h) { nav_h = 0; }
            $iframe = $(D.content);
            WIN_H = $(window).height();
            initIframe();

        };

    };
    
    Paper.sectionEdit = function(D) {
        getCache();
        var defaults = {

            "textWarp": "[tpl-section=warp]",
            "text": "[tpl-section=text]",
            "badge": "[tpl-section=badge]",
            "box": "[tpl-section=box]",
            "textarea": "[tpl-section=textarea]",
            "template": 'script[type="text/html"][template="section"]'
        };
        D = $.isPlainObject(D) ? $.extend(defaults, D) : defaults;
        var old_dom = null;
        return function() {

            old_dom = this;
            var text = "";
            var $this = $(this);
            var $warp = $this.closest(D.textWarp);
            var id = $warp.find(D.badge).data('id');
            var $text = $warp.find(D.text);
            var $box = $warp.find(D.box);
            $box.html($(D.template).html());
            var $textarea = $warp.find(D.textarea);
            cache.cache('id', id, function(index) {
                if (-1 === index) {
                    text = $text.text();
                } else {
                    var Obj = cache.get(index);
                    text = Obj.text;

                }
                $textarea.val($.trim(text));
                $textarea.text($textarea.val()).focus();
            });
        };

    };

    
    Paper.sectionSave = function(D) {
        getCache();
        var defaults = {
            'warp': '[tpl-section="warp"]',
            'box': '[tpl-section="box"]',
            'badge': '[tpl-section="badge"]',
            'textarea': '[tpl-section="textarea"]'
        };
        D = $.isPlainObject(D) ? $.extend(defaults, D) : defaults;
        var old_dom = null;
        return function() {

            var $warp = $(this).closest(D.warp);
            var text = $warp.find(D.textarea).val();
            text = $.trim(text);
            var id = $warp.find(D.badge).data('id');
            cache.cache('id', id, function(index) {
                if (-1 === index) {
                    cache.set({ 'text': text }, 'id', id);
                } else {
                    cache.update(index, { 'id': id, 'text': text });
                }
            });
            $warp.find(D.box).html('');
        };

    };

    
    Paper.textSectionEdit = function() {
        getCache();
        return function() {
            var text = "";
            var id = $(".section-id").text();
            var $textarea = $(".section-txt.section-txt-next");
            $textarea.text('');
            var arr = JSON.parse(section_edit_json);
            $.each(arr[id], function(index, value) {
                text += value.content;
            });
            cache.cache('id', id, function(index) {
                if (-1 === index) {
                    text = text;
                } else {
                    var Obj = cache.get(index);
                    text = Obj.text;
                }
                $textarea.val($.trim(text));
                $textarea.html($textarea.val());
            })
        }
    };
    
    Paper.textSectionSave = function() {
        getCache();
        return function() {
            var $this = $(this);
            var id = $this.parents('#layui-layer-section').find('.section-id').text();
            var text = $(".section-txt.section-txt-next").val();
            text = $.trim(text);
            cache.cache('id', id, function(index) {
                if (-1 === index) {
                    cache.set({ 'text': text }, 'id', id);
                } else {
                    cache.update(index, { 'id': id, 'text': text });
                }
            });
        }
    };

    
    Paper.tab = function() {
        var old_dom = null;
        return function(D) {
            if (old_dom === this) { return; }
            var defaults = {
                "li": "li",
                "section": '[tab="section"]',
                "active": 'active',
                "callback": function() {}
            };
            D = $.isPlainObject(D) ? $.extend(defaults, D) : defaults;
            old_dom = this;
            var $this = $(this);
            $this.parent().find(D.li).removeClass(D.active);
            $this.addClass(D.active);
            var id = $this.data('id');
            if (!id) { return; }
            var ids = id.toString().split(',');
            var $section = $(D.section);
            $section.hide();
            $section.each(function() {
                var $this = $(this);
                var id = $this.data('id');
                if ($.inArray(id.toString(), ids) !== -1) {
                    D.callback.call(this);
                    $this.show();
                }
            });

        };
    };


    
    Paper.toggle = function(show_fn, hide_fn) {
        var $section = $(this).closest('[eve-toggle=warp]').find('[eve-toggle=section]');
        if ($section.css("display") === "none") {
            if ($.isFunction(show_fn)) { show_fn(); }
            $section.show();
        } else {
            if ($.isFunction(hide_fn)) { hide_fn(); }
            $section.hide();
        }
    };

    System.Paper = Paper;
})(Report, jQuery);

(function(System, $) {
    'use strict';

})(Report, jQuery);


var flag = true;
var isAddLine = 'N'

function auxiliaryLine() {
    if (flag) {
        $('.addLine').html('取消辅助线');
        $('.red').css('border-bottom', '1px solid #f12828');
        $('.orange').css('border-bottom', '1px dashed #f39800');
        isAddLine = 'Y'
        return flag = false;
    } else {
        $('.addLine').html('使用辅助线');
        $('.red').css('border-bottom', '0px');
        $('.orange').css('border-bottom', '0px');
        isAddLine = 'N'
        return flag = true;
    }
}
$('.addLine').click(function() {
    auxiliaryLine();
})


$('.pagination > li > a').click(function(e) {
    e.preventDefault();
    location.href = $(this).attr('href') + '?isAddLine=' + isAddLine;
});


var theRequest = GetRequest()
if (theRequest.isAddLine === 'Y') {
    auxiliaryLine()
}

function GetRequest() {
    var url = location.search;
    var theRequest = new Object();
    if (url.indexOf("?") != -1) {
        var str = url.substr(1);
        strs = str.split("&");
        for (var i = 0; i < strs.length; i++) {
            theRequest[strs[i].split("=")[0]] = unescape(strs[i].split("=")[1]);
        }
    }
    return theRequest;
}