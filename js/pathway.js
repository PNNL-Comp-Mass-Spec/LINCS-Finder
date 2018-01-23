/**
 * 
 */
function KeggPathway(pid, sids, selector) {
    var _this = this;

    if(!d3) {
        console.log('d3 not installed');
        return;
    }
    console.log('pathway:' + pid);
    console.log('signature ids:' + sids);
    this.margin = {
        top:50,
        left:50,
        bottom:50,
        right:50,
        leftRight:100,
        topBottom:50
    }
    this.pathway_id = pid;
    this.sids = sids;
    this.sigData = {};
    this.sigGenes = {}

    
    ////////////////////////////////////////////////
    // manage the user input signature data
    
    var user_genes = localStorage.getItem('user_genes')
    if (user_genes != null) {
        // kegg id <--> gene symbol for user data
        this.kegg_genes = JSON.parse(localStorage.getItem('kegg_genes'));
        this.userGeneSym2KeggId = {};
        this.geneIdToSymbol ={};
        this.kegg_genes.forEach(function(d){
            var kid = d.gene.split(':')[1];
            _this.userGeneSym2KeggId[d.symbol] = kid;
            _this.geneIdToSymbol[kid] = d.symbol;
        });

        this.sigData['User_Input_Sig'] = {cid:"User_Input_Sig", sig_id:"User_Input_Sig", pert_iname:"User_Input_Sig"}
        var userGenes = JSON.parse(user_genes);
    
        this.sigGenes['User_Input_Sig'] = {}
        this.sigGenes['User_Input_Sig']['dn'] = []
        this.sigGenes['User_Input_Sig']['up'] = []
        this.sigGenes['User_Input_Sig']['dn_ezids'] = []
        this.sigGenes['User_Input_Sig']['up_ezids'] = []
     
        userGenes.dn.forEach(function(sym){
            _this.sigGenes['User_Input_Sig'].dn.push(sym);
            _this.sigGenes['User_Input_Sig'].dn_ezids.push(_this.userGeneSym2KeggId[sym]);
        });
        userGenes.up.forEach(function(sym){
            _this.sigGenes['User_Input_Sig'].up.push(sym);
            _this.sigGenes['User_Input_Sig'].up_ezids.push(_this.userGeneSym2KeggId[sym]);
        });
        console.log(this.kegg_genes)
        console.log(this.sigGenes)
    }
    
    ////////////////////////////////////////////////
    
    this.svgWidth = window.innerWidth-this.margin.leftRight;
    this.svgHeight = window.innerHeight-this.margin.topBottom;
    this.selector = selector.append('svg');
    
    this.setup();
}

function extend(obj, src) {
    Object.keys(src).forEach(function(key) { obj[key] = src[key]; });
    return obj;
}

KeggPathway.prototype = {

    getSigData: function(sig_info) {
        var _this = this;
        sig_info.forEach(function(d){
            _this.sigData[d.cid] = d;
        });
    },

    getGeneData: function(gene_info) {
        this.sigGenes = extend(gene_info, this.sigGenes);
    },

    infoButtonEvt: function() {
        var content = '<table class="ui table">';
        var body = [];
        var sids = Object.keys(this.sigData);
        var _this = this;
        content += '<thead><tr><th>Attribute</th>';
        sids.forEach(function(d){
            content += '<th>'+d+'</th>'
        });
        content += '</tr></thead>';
        var attKeys = Object.keys(_this.sigData[sids[sids.length-1]]);
        content += '<tbody>'
        for(var i = 0; i<attKeys.length; i++){
            content+='<tr><td>'+attKeys[i]+'</td>';
            for(var j =0; j<sids.length; j++){
                content+='<td>'+_this.sigData[sids[j]][attKeys[i]]+'</td>';
            }
            content+='</tr>'
        }
        content += '</tbody></table>';
        
        $('#nodeInfo').popup({
            on:'click',
            position:'bottom left',
            html:content
        });
    },
        setup : function() {
            var _this = this;
            var fetchGenePromise = new Promise(function(resolve, reject) { 
             _this.fetchGenesFromNodes(_this.sids, resolve, reject);
            });
            var fetchKgmlPromise = new Promise(function(resolve, reject) { 
                _this.fetchKGML(resolve, reject);
            });
            
            Promise.all([fetchKgmlPromise, fetchGenePromise]).then(values => { 
                _this.initPathway(values[0]);
                _this.getGeneData(values[1].genes);
                _this.setColors();
                _this.getSigData(JSON.parse(values[1].sigs));
                _this.infoButtonEvt();
            });
        },
            
        fetchKGML : function(resolve, reject) {
            var _this = this;
            var kgmlUrl = "https://adbio.pnnl.gov/bioviz/services/kgml/" + _this.pathway_id;
            _this.ajaxCall(kgmlUrl,'POST',{},
                    function(data, status) {
                resolve(data);
            }, function(data, status, c) {
                reject(Error("Network Error"));
            });
        },
        
        fetchGenesFromNodes : function(sids, resolve, reject) {
            var lincsGenes = {};
            // var nodesForApi = [];
            // sids.forEach(function(sid) {
            //     if (localStorage.getItem(node)) {
            //         lincsGenes[node] = localStorage.getItem(node);
            //     } else {
            //         nodesForApi.push(node);
            //     }
            // });
            
            // if (nodesForApi.length > 0) {
            this.fetchGenesFromApi(sids, resolve);
            // } else {
            //     resolve(lincsGenes);
            // }
        },
        
        // fetchGenesFromApi : function(nodes, lincsGenes, resolve) {
        //   var _this = this;
        //     $.ajax({
        //         url:'/lincs/genes/nodes/'+(nodes.join && nodes.join(',')||nodes)
        //     }).done(function(data){
        //         var nodeGenes = JSON.parse(data);
        //         // console.log(nodeGenes);
        //         nodeGenes.nodes.forEach(function(d){
        //             lincsGenes[d._id['$oid']] = d;
        //         });
        //         nodeGenes.gene_info.filter(function(d){return d.gene_info[0].kegg;}).forEach(function(d){
        //           var tmp = d.gene_info[0].kegg.split(':');
        //           _this.geneIdToInfo[tmp[1]]=d.gene_info[0];
        //         });
        //         resolve(lincsGenes);
        //     });
        // },

        fetchGenesFromApi : function(sig_id, resolve) {
            var _this = this;
            $.ajax({
                url:'/lincs/signature?sig_id='+sig_id[0],
                method: 'GET',
                dataType: "jsonp",
                jsonpCallback: "logResults"
            }).done(function(data) {
                // console.log(JSON.parse(data.sig)[0])
                resolve(data);
                // fillTable('pos', JSON.parse(data.pos))
                // fillTable('neg', JSON.parse(data.neg))
                // $('table#pos').DataTable();
                // $('table#neg').DataTable();
            });
            // $.ajax({
            //     url:'/lincs/genes/nodes/'+(sig_ids.join && sig_ids.join(',')||sig_ids)
            // }).done(function(data){
            //     var nodeGenes = JSON.parse(data);
            //     // console.log(nodeGenes);
            //     nodeGenes.nodes.forEach(function(d){
            //         lincsGenes[d._id['$oid']] = d;
            //     });
            //     nodeGenes.gene_info.filter(function(d){return d.gene_info[0].kegg;}).forEach(function(d){
            //       var tmp = d.gene_info[0].kegg.split(':');
            //       _this.geneIdToInfo[tmp[1]]=d.gene_info[0];
            //     });
            //     resolve(lincsGenes);
            // });
        },

        convertXml : function(pathway) {
            var _this = this;
            var graphics = function(selector) {
                var obj = {
                    bgcolor : selector.attr('bgcolor'),
                    fgcolor : selector.attr('fgcolor'),
                    height : selector.attr('height'),
                    name : selector.attr('name'),
                    type : selector.attr('type'),
                    width : selector.attr('width'),
                    x : selector.attr('x'),
                    y : selector.attr('y'),
                    coords : selector.attr('coords')
                };
                return obj;
            };
            
            this.pathwayEntry = [];
            this.geneName2entry = {};
            var sids = Object.keys(this.sigData);
            var scales = d3.scaleLinear().domain([0, sids.length+1]);
            
            var entry = function(selector) {
                var obj = {
                        id : selector.attr('id'),
                        link : selector.attr('link'),
                        name : selector.attr('name'),
                        type : selector.attr('type'),
                        graphics : graphics(selector.find('graphics')),
                        dn : [],
                        up : [],
                        color: function(d) {
                            if(this.dn.length*this.up.length > 0)
                                return 'green';
                            else if (this.dn.length > 0) {
                                var c = scales.range(['white', 'blue']);
                                return c(this.dn.length);
                            } else if (this.up.length > 0){
                                var c = scales.range(['white', 'red']);
                                return c(this.up.length);
                            }
                            else return 'white';
                        }
                };
                id2entry = {};
                id2entry[obj.id] = obj;
                _this.pathwayEntry.push(id2entry);
                
                if (obj.name.indexOf('hsa:') >= 0) {
                    obj.name.split(' ').forEach(function(n){
                        if (n.indexOf('hsa:') >= 0) {
                            if (!_this.geneName2entry.hasOwnProperty(n))
                                _this.geneName2entry[n] = [];
                            _this.geneName2entry[n].push(obj);
                        }
                    });
                    
                }
                
                return obj;
            };
            
            this.pathwayJson = pathway.map(function(d) {
                var obj = {
                        image : $(this).attr('image'),
                        link : $(this).attr('link'),
                        name : $(this).attr('name'),
                        number : $(this).attr('number'),
                        org : $(this).attr('org'),
                        title : $(this).attr('title'),
                        entry : $(this).find('entry').map(function(d) {
                            return entry($(this));
                        })
                };
                return obj;
            })[0];
        },
        ajaxCall : function(url,type,data,success,fail){
            var ajaxObj = {
                    url:url,
                    type:type||'GET'
            };
            if(data)
                ajaxObj['data'] = data;

            $.ajax(ajaxObj).done(success).fail(fail);
        },
        onAjaxError:function(data, status, c) {
            console.log(data);
            console.log(status);
            console.log(c);
        },
        initPathway : function(data) {
            var _this = this;
            //store kgml to div
            $('#kegg').html(data);
            var pathway = $('#kegg').find('pathway');

            //convert kgml to json object
            _this.convertXml(pathway);
            _this.initView();
            _this.refresh();
        },
        setColors : function() {
            var _this = this;
            Object.keys(this.sigGenes).forEach(function(id){
                var node = _this.sigGenes[id];
                node.dn_ezids.forEach(function(d, i) {
                    var gene_id = ''+d;
                    var gene_sym = node.dn[i];
                    var entries = _this.geneName2entry['hsa:'+d];
                    if (entries) {
                        entries.forEach(function(e){
                            e.dn.push({nodeId:id, geneId:gene_id, geneSym:gene_sym});
                            _this.geneIdToSymbol[gene_id] = gene_sym
//                          e.color = 'blue';
//                          if(e.dn.length*e.up.length > 0) e.color = 'green'; 
                        })
                    }
                });
                node.up_ezids.forEach(function(d, i) {
                    var gene_id = ''+d;
                    var gene_sym = node.up[i];
                    var entries = _this.geneName2entry['hsa:'+d];
                    if (entries) {
                        entries.forEach(function(e){
                            e.up.push({nodeId:id, geneId:gene_id, geneSym:gene_sym});
                            _this.geneIdToSymbol[gene_id] = gene_sym
//                          e.color = 'red';
//                          if(e.dn.length*e.up.length > 0) e.color = 'green'; 
                        })
                    }
                });
            });
            this.refresh();
        },
        initView:function(){
            var _this = this;
            //zoom and pan function
            function zoomed() {
                if(d3.event.transform){
                    _this.selector.select('g').attr('transform', d3.event.transform);
                }
            }
            _this.zoom = d3.zoom().scaleExtent([ -10, 10 ]).on("zoom", zoomed);
            
            var tipTimeout = {};
            var tip = function(d, _target, _class, funct) {
                var a = {
                    show: function() {
                        if (tipTimeout[d.id]) {
                            clearTimeout(tipTimeout[d.id]);
                            tipTimeout[d.id] = undefined;
                            return;
                        }
                        var nodeIds = d.dn.map(function(f){return f.nodeId;})
                            .concat(d.up.map(function(f){return f.nodeId;}))
                            .filter(function(item, pos,a) { return a.indexOf(item) == pos;});
                        var upGenes = d.up;
                        var dnGenes = d.dn;
                        var nodeObj=[];
                        // var columns = ['node_id'];
                        var columns = [];
                        nodeIds.forEach(function(f){
                            var tmp = {
                                node_id:f
                            };
                            upGenes.forEach(function(g){
                                if(g.nodeId == f){
                                    tmp[g.geneId] = 'up';
                                    columns.push(g.geneId);
                                }
                            });
                            dnGenes.forEach(function(g){
                                if(g.nodeId == f){
                                    if(tmp.hasOwnProperty(g.geneId)) tmp[g.geneId] = 'warning';
                                    else tmp[g.geneId] = 'down';
                                    columns.push(g.geneId);
                                }
                            })
                            nodeObj.push(tmp);
                        });
                        var heads = d.name.split(' ').map(function(g){return g.split(':')[1];})
                                .filter(function(f){return columns.indexOf(f)>=0;});
                        var html = '<table class="ui striped celled table"><thead><tr><th>Signatures</th>';//<th class="center aligned">'
                        heads.forEach(function(f){
                            html+='<th class="center aligned">'+_this.geneIdToSymbol[f]+'</th>';
                            // html+='<th class="center aligned">'+d.name+'('+f+')'+'</th>';
                        });
                        // +heads.join('</th><th class="center aligned">')+'</th></tr></thead>'
                        html += '</tr></thead><tbody>';
                        nodeObj.forEach(function(f){
                            var popupContent='',
                                tdataContent=f.node_id;
                            if(_this.sigData[f.node_id]){
                                popupContent = _this.sigData[f.node_id].sig_id;
                                tdataContent = _this.sigData[f.node_id].cid;
                            }
                            html += '<tr class="center aligned" data-content="' + popupContent + '"><td>'
                                + tdataContent + '</td>';
                            heads.forEach(function(g){
                                html += '<td><i class="large '+(f[g]==='up'?'outline red':f[g]==='down'?'outline blue':f[g]==='warning'?'green':'outline')+' arrow circle '+f[g]+' icon"></i></td>';
                            });
                            html += '</tr>';
                        });
                        html += '</tbody></table>';
                        var position = {
                            top:function(){
                                var targetTop = _target.getBoundingClientRect().top,
                                    targetHeight = _target.getBoundingClientRect().height,
                                    thisTop = this.getBoundingClientRect().top,
                                    thisHeight = this.getBoundingClientRect().height,
                                    winHeight = window.innerHeight || window.clientHeight;
                                var pos = targetTop;
                                if(targetTop+thisHeight >= winHeight)
                                    pos -= (thisHeight + targetHeight);
                                if(pos < 0) pos = 0;
                                return (pos)+'px';
                            },
                            left:function(){
                                var targetLeft = _target.getBoundingClientRect().left,
                                    targetWidth = _target.getBoundingClientRect().width,
                                    thisLeft = this.getBoundingClientRect().left,
                                    thisWidth = this.getBoundingClientRect().width,
                                    winWidth = window.innerWidth || window.clientWidth;
                                var pos = targetLeft+targetWidth;
                                if(pos + thisWidth >= winWidth)
                                    pos -= (thisWidth + targetWidth);
                                if(pos < 0) pos = 0;
                                return (pos)+'px';
                            }
                        }
                        
                        if(heads.length > 0){
                            var tipContainer = d3.select('body').append('div').attr('class', 'ui segment tip')
                                    .attr('id','tip_'+d.id)
                                    .style('position','absolute')
                                    .html(html).on('mouseover',function(){
                                        if(tipTimeout[d.id]){
                                            clearTimeout(tipTimeout[d.id]);
                                            tipTimeout[d.id] = undefined;
                                        }
                                    }).on('mouseout',function(){
                                        var _evt = this;
                                        tipTimeout[d.id] = setTimeout(function(){
                                            d3.select(_evt).transition().duration(500).style('transform','scale(0)').remove();
                                            if(tipTimeout[d.id]){
                                                clearTimeout(tipTimeout[d.id]);
                                                tipTimeout[d.id] = undefined;
                                            }
                                        },500);
                                    });
                            
                            tipContainer.styles({
                                'top':position.top,
                                'left':position.left,
                                'max-height':window.innerHeight+'px',
                                'overflow':'auto'
                            });
                            $('#tip_'+d.id).find('tr').popup();
                        }
                    },
                    hide:function(){
                        tipTimeout[d.id] = setTimeout(function(){
                            d3.select('body').selectAll('#tip_'+d.id).transition().duration(500).style('transform','scale(0)').remove();
                            if(tipTimeout[d.id]){
                                clearTimeout(tipTimeout[d.id]);
                                tipTimeout[d.id] = undefined;
                            }
                        },500);
                    }
                };
                a[funct]();
            }

            var entriesRefresh = function(d) {
                var obj = {};
                var key = 'fill';
                if (['line','poly'].indexOf(d.graphics.type)/*.indexOf('line')*/ >= 0) {
                    key = 'stroke';
                }
//              obj[key] = Array.isArray(d.color)?d.color[0]:d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:'white';
                obj[key] = d.color();

                d.dn.forEach(function(g) {
                    if (g.nodeId=='User_Input_Sig'){
                        // if (!('stroke' in obj)) {
                        if (Object.keys(obj).indexOf('stroke')<0) {
                            obj['stroke'] = 'red'
                        }
                        obj['stroke-width'] = 5
                    }
                });

                d.up.forEach(function(g) {
                    if (g.nodeId=='User_Input_Sig'){
                        // if (!('stroke' in obj)) {
                        if (Object.keys(obj).indexOf('stroke')<0) {
                            obj['stroke'] = 'red'
                        }
                        obj['stroke-width'] = 5
                    }
                });
                
                return obj;
            };
            var entriesAppend = function(d) {
                if (d.graphics.type === 'rectangle' || d.graphics.type === 'roundrectangle')
                    return document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                if (d.graphics.type === 'circle')
                    return document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                if ((d.graphics.type === 'line' || d.graphics.type === 'poly'))// && d3.keys(d.values).length > 0)
                    return document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
                //return div element if the other types not found
                return document.createElement('div');
            };
            var entriesAttributes = function(d) {
                var type = d.graphics.type;
                //attributes for rectangle
                if (type.indexOf('rect') >= 0) {
                    var obj = {
                            id : d.id,
                            x : Math.ceil(d.graphics.x - d.graphics.width / 2) + 1,
                            y : Math.floor(d.graphics.y - d.graphics.height / 2) + 1,
                            width : (+d.graphics.width - 1),
                            height : (+d.graphics.height - 1),
//                          fill : d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:undefined
//                          || d.color[0],
                            opacity : .5
                    };
                    //add rounded corners
                    if (type.indexOf('round') >= 0 && _this.pathwayJson.title != d.graphics.name) {
                        obj['rx'] = 8;
                        obj['ry'] = 8;
                    }
                    return obj;
                } else if (type.indexOf('circle') >= 0) {
                    //attributes if circle
                    return {
                        id : d.id,
                        cx : +d.graphics.x + .5,
                        cy : +d.graphics.y + .5,
                        r : (+d.graphics.width / 2),
//                      fill : d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:undefined
//                      || d.color[0],
                        opacity : .5
                    };
                } else if (type.indexOf('line') >= 0) {
                    //attributes if line
                    return {
                        id : d.id,
                        points:d.graphics.coords,
                        'stroke-linejoin':'round',
                        'stroke-linecap':'round',
                        fill : 'none',
//                      stroke : d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:undefined
//                      || d.color[0],
                        'stroke-width' : 8,
                        opacity : .5
                    };
                } else if (type.indexOf('poly') >= 0) {
                    //attributes if line
                    return {
                        id : d.id,
                        points:d.graphics.coords,
                        'stroke-linejoin':'round',
                        'stroke-linecap':'round',
                        fill : 'none',
//                      stroke : d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:undefined
//                      || d.color[0],
                        'stroke-width' : 8,
                        opacity : .5
                    };
                }
            };

            var svg = this.selector.attrs({
                'height' : _this.svgHeight,
                'width' : _this.svgWidth,
            });
            
            svg.call(_this.zoom);
            svg.append('filter').attr('id','grayscale').append('feColorMatrix').attr('type','saturate').attr('values','0');
            var g = svg.append('g');
            g.append('image').attrs({
                href : _this.pathwayJson.image,
                x : 0,
                y : 0,
                //filter:"url(#grayscale)"
            });

            _this.refresh = function() {
                // refreshing the view with data
                var entries = g.selectAll('.entry').data(_this.pathwayJson.entry)
                .attrs(entriesRefresh);

                // creating the dom elelments
                var groupEntry = entries.enter().append('g');

                groupEntry.append(entriesAppend)
                .classed('entry', true)
                .style('cursor', 'pointer')
                .attrs(entriesAttributes)
                .on('click', function(d) {
                    _this.generateModal($('body'), d);
                }).on('mouseover', function(d) {
                    d3.select(this).attr('opacity', 1);
                    tip(d,this,_this,'show');
                }).on('mouseout', function(d) {
                    d3.select(this).attr('opacity', .5);
                    
                      tip(d,this,_this,'hide');
                });
//              };
                groupEntry.append(function(d){
                    if (d.graphics.type === 'line' && d3.keys(d.values).length > 0)
                        return document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
                    //return div element if the other types not found
                    return document.createElement('div');
                }).attrs(function(d){
                    if (d.graphics.type.indexOf('line') >= 0) {
                        //attributes if line
                        return {
                            id : d.id,
                            points:d.graphics.coords,
                            'stroke-linejoin':'round',
                            'stroke-linecap':'round',
                            fill : 'none',
//                          stroke : d.color[_this.timeIndex]?d.color[_this.timeIndex][_this.groupIndex]:undefined
//                          || d.color[0],
                            'stroke-width' : 8,
                            opacity : .5
                        };
                    }
                })
            }
        }
}
