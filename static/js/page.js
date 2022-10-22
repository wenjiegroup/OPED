    $(function () {
        // ajax分页
		// 打开页面加载第一页
        let now_page = 1;
		page_click()
        $('.first_page button').removeClass('layui-btn-primary').addClass('layui-btn-disabled');
        $('.now_page button').first().removeClass('layui-btn-primary').addClass('page_this');
        //上一页
        $('.first_page').click(function () {
            now_page -= 1;
            if (now_page < 1) {
                now_page = 1;
                return false
            } else {
                $('.page_this').parent().prev().click();
            }
        });
        //下一页
        $('.last_page').click(function () {
            let num_pages = $('.now_page button').last().text();
            now_page += 1;
            if (now_page > parseInt(num_pages)) {
                now_page -= 1;
                return false
            } else {
                $('.page_this').parent().next().click();
            }
        });
        //切换页
        $('.now_page').click(function () {
            now_page = parseInt($(this).children('button').text());
            $('.now_page button').removeClass('page_this').addClass('layui-btn-primary');
            $(this).addClass('page_this');
            $(this).children('button').removeClass('layui-btn-primary').addClass('page_this');
            page_click()
        });

        function page_click() {
            let page_form = $('#page');
            $.ajax({
                type: 'get',
                url: page_form.attr('action'),
                data: {page: now_page},
                success: function (data) {
                    $('#tbody tr').remove();
                    $('#num_pages').html('共' + data.num_pages + '页');
                    if (data.has_previous === true) {
                        $('.first_page button').removeClass('layui-btn-disabled').addClass('layui-btn-primary');
                    } else {
                        $('.first_page button').removeClass('layui-btn-primary').addClass('layui-btn-disabled')
                    }
                    if (data.has_next === true) {
                        $('.last_page button').removeClass('layui-btn-disabled').addClass('layui-btn-primary');
                    } else {
                        $('.last_page button').removeClass('layui-btn-primary').addClass('layui-btn-disabled');
                    }
                    $.each(data, function (index, user) {
                        let a = '<td>';
                        let b = '</td>';
                        let body = a + peg.username + b + a + peg.is_rank + b + a + user.sex + b + a + user.phone + b + a + user.email + b + a + user.last_login + b + a + user.date_joined + b + a + user.update_time + b;
                        $('#tbody').append('<tr>' + body + '</tr>');
                    });
                }
            })
        }
    })